#include "vbscf/optimization/backends/ciah.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "vbscf/derivatives/hessian/coupled/structure.hpp"
#include "vbscf/derivatives/hessian/exact/operator.hpp"
#include "vbscf/derivatives/hessian/context/accepted_point.hpp"
#include "vbscf/diagnostics/performance/neo_trace.hpp"
#include "vbscf/optimization/driver/options.hpp"
#include "vbscf/optimization/driver/result.hpp"
#include "vbscf/optimization/globalization/line_search.hpp"
#include "vbscf/optimization/neo/ciah.hpp"
#include "vbscf/optimization/neo/globalization.hpp"
#include "vbscf/optimization/trust_region/retraction.hpp"

namespace xmvb::vb::optimizer_detail {
namespace {

// Initial policies follow PySCF newton_casscf.update_orb_ci. They control
// approximate microsteps, not the VBSCF outer accuracy contract.
constexpr double max_component = 0.03;
constexpr double ah_start_tolerance = 500.0;
constexpr double gradient_trust = 3.0;
constexpr int ah_start_cycle = 3;
constexpr int keyframe_interval = 5;
constexpr int base_micro_cycles = 10;

struct Keyframe {
  Eigen::VectorXd gradient;
  double energy = 0.0;
};

}  // namespace

bool build_ciah_macro_step(
    VbScfObjective* objective,
    const SparseParameterLayout& layout,
    const VbScfOptimizerOptions& options,
    const OrbitalChart& chart,
    const NonredundantRetractionMetric& metric,
    ExactHvpOperator& hessian,
    const StructureTangentOperator& structure,
    const Eigen::VectorXd& gradient,
    const Eigen::VectorXd& orbital_guess_packed,
    const Eigen::MatrixXd& structure_guess,
    bool use_structure_response,
    double energy,
    double* trust_radius,
    NeoIterationRecord* record,
    VbScfOptimizerResult* result,
    CiahMacroStep* accepted) {
  const auto context = objective->second_order_context();
  if (context->selected_state_indices.size() != 1) {
    throw std::invalid_argument("coupled CIAH requires a single selected state");
  }
  const auto orbitals = objective->input().orbital_preparation_input;
  const Eigen::Index no = chart.reduced_size();
  const Eigen::Index ns = use_structure_response ? structure.tangent_size() : 0;
  const Eigen::Index dimension = no + ns;
  record->model_dimension = static_cast<int>(dimension);
  record->source_gradient_norm = gradient.norm();
  record->source_physical_gradient_norm = metric.dual_norm(gradient);

  const auto dual_norm = [&](const Eigen::VectorXd& g) {
    return std::hypot(metric.dual_norm(g.head(no)),
        ns ? structure.diagonal_dual_norm(g.tail(ns)) : 0.0);
  };
  const auto action = [&](const Eigen::VectorXd& x) {
    Eigen::VectorXd image(dimension);
    if (ns == 0) {
      image = hessian.apply_orbital(x);
    } else {
      const auto orbital = hessian.apply_orbital_coupling(x.head(no));
      const auto ci = structure.apply_coupling_coordinates(x.tail(ns));
      image.head(no) = orbital.orbital_hessian +
          hessian.apply_structure_coupling_adjoint(
              ci.coefficient_response, ci.adjoint_multipliers);
      image.tail(ns) = structure.project_coordinates(
          orbital.scaled_structure_forcing) + ci.hessian_coordinates;
      ++record->structure_response_actions;
    }
    ++record->orbital_hvp_actions;
    ++record->coupled_block_actions;
    return image;
  };
  const auto metric_action = [&](const Eigen::VectorXd& x) {
    Eigen::VectorXd image(dimension);
    image.head(no) = metric.apply(x.head(no));
    if (ns) image.tail(ns) = structure.apply_metric_coordinates(x.tail(ns));
    return image;
  };
  const auto precondition = [&](const Eigen::VectorXd& g, double eigenvalue) {
    Eigen::VectorXd direction(dimension);
    if (options.orbital_preconditioner == OrbitalPreconditioner::Identity) {
      direction.head(no) = g.head(no);
    } else if (options.orbital_preconditioner == OrbitalPreconditioner::OneElectron) {
      // The one-electron model is a positive block preconditioner, not a
      // signed approximation to the complete AH operator. Restore its
      // unshifted inverse; the complete physical metric stays in metric_action.
      direction.head(no) = chart.apply_inverse_reduced_block_preconditioner(g.head(no));
    } else {
      direction.head(no) = chart.apply_inverse_augmented_hessian_diagonal(
          g.head(no), eigenvalue - 1.0e-8);
    }
    if (ns) direction.tail(ns) = structure.apply_inverse_augmented_hessian_diagonal(
        g.tail(ns), eigenvalue - 1.0e-8);
    return direction;
  };
  Eigen::VectorXd initial_gradient = Eigen::VectorXd::Zero(dimension);
  initial_gradient.head(no) = gradient;
  NeoProblem problem(initial_gradient, action, metric_action, dual_norm);
  CiahOptions ciah_options;
  ciah_options.residual_preconditioner = precondition;
  double energy_scale = 0.0;
  for (const double state_energy : context->selected_state_energies)
    energy_scale = std::max(energy_scale, std::abs(state_energy));
  ciah_options.operator_relative_accuracy =
      context->structure_solve_accuracy.response_backward_error_tolerance(energy_scale);
  Eigen::VectorXd seed = initial_gradient;
  if (orbital_guess_packed.size() == static_cast<Eigen::Index>(layout.size())) {
    seed.head(no) = chart.project_vector(orbital_guess_packed).reduced_gradient;
  }
  if (ns && structure_guess.rows() == context->selected_state_eigenvectors.rows() &&
      structure_guess.cols() == 1) {
    seed.tail(ns) = structure.project_coordinates(
        std::sqrt(2.0 * context->normalized_state_weights[0]) * structure_guess);
  }
  CiahWorkspace workspace(problem, seed, ciah_options);

  const Eigen::MatrixXd c0 = context->selected_state_eigenvectors;
  const double coordinate_scale = std::sqrt(2.0 * context->normalized_state_weights[0]);
  const auto coefficients = [&](const Eigen::VectorXd& displacement) {
    Eigen::MatrixXd c = c0;
    if (ns) c += structure.coefficient_response_coordinates(displacement.tail(ns));
    return c;
  };
  const auto evaluate_keyframe = [&](const Eigen::VectorXd& displacement) {
    const auto parameters = build_nonredundant_lifted_trial_parameters(
        orbitals, chart, layout, displacement.head(no));
    const Eigen::MatrixXd c = coefficients(displacement);
    const auto exact = objective->evaluate_coupled_keyframe(parameters, c);
    Keyframe keyframe;
    keyframe.gradient.resize(dimension);
    keyframe.gradient.head(no) = chart.project_reduced_gradient(exact.orbital_gradient);
    if (ns) {
      // C_raw/norm is the finite overlap-normalized coefficient retraction.
      // Its derivative contributes 1/norm; simply projecting (H-E M)C
      // would compare gradients from different coefficient coordinates.
      const auto normalized = exact.normalized_structure_coefficients.col(0);
      const double norm = c.col(0).dot(normalized) / normalized.squaredNorm();
      if (!(norm > 0.0) || !std::isfinite(norm)) {
        throw std::runtime_error("CIAH keyframe coefficient normalization is singular");
      }
      keyframe.gradient.tail(ns) = structure.project_coordinates(
          (coordinate_scale / norm) * exact.structure_residuals);
    }
    keyframe.energy = exact.energy;
    ++record->keyframes;
    return keyframe;
  };

  Eigen::VectorXd displacement = Eigen::VectorXd::Zero(dimension);
  Eigen::VectorXd hessian_displacement = displacement;
  Eigen::VectorXd predicted_gradient = initial_gradient;
  Eigen::VectorXd last_step = displacement;
  Eigen::VectorXd descent_displacement, descent_image, descent_last_step;
  double best_predicted_reduction = 0.0;
  // Keyframe corrections need not preserve descent from the macro origin.
  // Retain only prefixes left in place by the existing microstep checks.
  // A negative origin slope allows nonlinear backtracking to find descent;
  // no extra Hessian action is needed to retain its frozen model image.
  const auto retain_descent_prefix = [&] {
    const double slope = initial_gradient.dot(displacement);
    const double reduction = -slope - 0.5 * displacement.dot(hessian_displacement);
    if (slope < 0.0 && reduction > best_predicted_reduction) {
      best_predicted_reduction = reduction;
      descent_displacement = displacement;
      descent_image = hessian_displacement;
      descent_last_step = last_step;
    }
  };
  double keyframe_norm = dual_norm(initial_gradient);
  // newton_casscf.kernel passes 0.3 * its outer gradient tolerance to
  // update_orb_ci; the microstep/keyframe tests use that inner tolerance.
  const double micro_gradient_tolerance = 0.3 * options.gradient_tolerance;
  int micro_steps = 0, since_keyframe = 0;
  CiahResult candidate;
  while (true) {
    retain_descent_prefix();
    candidate = workspace.next(predicted_gradient);
    ++record->micro_iterations;
    const bool ready = candidate.converged || candidate.exhaustion ||
        (candidate.basis_size >= ah_start_cycle &&
         candidate.residual_norm < ah_start_tolerance);
    if (!ready) continue;

    Eigen::VectorXd step = candidate.step;
    Eigen::VectorXd image = candidate.hessian_step;
    double largest = chart.maximum_rotation_component(step.head(no));
    if (ns) largest = std::max(largest,
        structure.maximum_coefficient_component_coordinates(step.tail(ns)));
    if (largest > max_component) {
      const double scale = max_component / largest;
      step *= scale;
      image *= scale;
    }
    const Eigen::VectorXd previous_gradient = predicted_gradient;
    const Eigen::VectorXd previous_step = last_step;
    displacement += step;
    hessian_displacement += image;
    predicted_gradient += image;
    last_step = step;
    ++micro_steps;
    ++since_keyframe;
    const double predicted_norm = dual_norm(predicted_gradient);
    diagnostics::neo_trace("ciah_micro",
        "basis=%d,micro=%d,omega=%.17g,alpha=%.17g,residual=%.17g,gradient=%.17g",
        candidate.basis_size, micro_steps, candidate.eigenvalue,
        candidate.reference_weight, candidate.residual_norm, predicted_norm);
    if (micro_steps > 3 && predicted_norm > gradient_trust * keyframe_norm) {
      displacement -= step;
      hessian_displacement -= image;
      predicted_gradient = previous_gradient;
      last_step = previous_step;
      ++record->rejected_trial_count;
      break;
    }
    const int micro_limit = std::max(base_micro_cycles,
        base_micro_cycles - static_cast<int>(2.0 * std::log(keyframe_norm + 1.0e-7)));
    if (micro_steps >= micro_limit ||
        predicted_norm < 0.3 * micro_gradient_tolerance || candidate.exhaustion) break;

    const double length = std::sqrt(std::max(0.0,
        displacement.dot(metric_action(displacement))));
    const double interval = std::max<double>(keyframe_interval,
        keyframe_interval - std::log(length + 1.0e-7));
    if (since_keyframe >= interval || predicted_norm < keyframe_norm / gradient_trust) {
      const Keyframe exact = evaluate_keyframe(displacement);
      const double actual_norm = dual_norm(exact.gradient);
      const double mismatch = dual_norm(exact.gradient - predicted_gradient);
      diagnostics::neo_trace("ciah_keyframe",
          "actual=%.17g,predicted=%.17g,mismatch=%.17g,energy=%.17g",
          actual_norm, predicted_norm, mismatch, exact.energy);
      if (mismatch < gradient_trust * predicted_norm ||
          actual_norm < gradient_trust * micro_gradient_tolerance) {
        predicted_gradient = exact.gradient;
        keyframe_norm = actual_norm;
        since_keyframe = 0;
      } else {
        displacement -= step;
        hessian_displacement -= image;
        predicted_gradient = previous_gradient;
        last_step = previous_step;
        ++record->rejected_trial_count;
        break;
      }
    }
  }

  retain_descent_prefix();
  const double final_slope = initial_gradient.dot(displacement);
  const double final_reduction =
      -final_slope - 0.5 * displacement.dot(hessian_displacement);
  bool restored_prefix = false;
  // Preserve legitimate finite negative-curvature steps with positive model
  // reduction, even if their initial slope is uphill.
  if (final_reduction <= 0.0 &&
      best_predicted_reduction > 0.0) {
    diagnostics::neo_trace("ciah_descent_prefix",
        "discarded_linear=%.17g,discarded_prediction=%.17g,prefix_prediction=%.17g",
        final_slope, final_reduction, best_predicted_reduction);
    displacement = std::move(descent_displacement);
    hessian_displacement = std::move(descent_image);
    last_step = std::move(descent_last_step);
    restored_prefix = true;
  }

  accepted->parameters = build_nonredundant_lifted_trial_parameters(
      orbitals, chart, layout, displacement.head(no));
  if ((accepted->parameters - layout.pack(orbitals)).norm() <=
      std::numeric_limits<double>::epsilon() * layout.pack(orbitals).norm()) {
    result->termination_reason = "ciah_zero_macro_step";
    return false;
  }
  // The frozen CIAH microsteps do not replace the VBSCF energy acceptance
  // contract. Reuse their accumulated H image to screen a shorter macro
  // displacement, without repeating any HVP or building another Hessian.
  double scale = 1.0;
  const double linear = initial_gradient.dot(displacement);
  const double quadratic = displacement.dot(hessian_displacement);
  while (true) {
    accepted->trial = objective->evaluate_trial_with_structure_guess(
        accepted->parameters, coefficients(scale * displacement));
    const auto globalization = globalize_neo_trial(energy, accepted->trial.energy,
        scale * linear, scale * scale * quadratic, false,
        *trust_radius, options.minimum_step_size);
    diagnostics::neo_trace("ciah_macro_trial",
        "scale=%.17g,linear=%.17g,quadratic=%.17g,current_energy=%.17g,trial_energy=%.17g,predicted=%.17g,actual=%.17g,rho=%.17g,accepted=%d",
        scale, scale * linear, scale * scale * quadratic, energy,
        accepted->trial.energy, globalization.predicted_reduction,
        globalization.actual_reduction, globalization.rho,
        static_cast<int>(globalization.accepted));
    if (globalization.accepted) {
      *trust_radius = globalization.next_radius;
      break;
    }
    ++record->rejected_trial_count;
    scale *= 0.5;
    if (scale * metric.norm(displacement.head(no)) <= options.minimum_step_size) {
      result->termination_reason = "ciah_macro_energy_rejected";
      return false;
    }
    accepted->parameters = build_nonredundant_lifted_trial_parameters(
        orbitals, chart, layout, scale * displacement.head(no));
  }
  displacement *= scale;
  hessian_displacement *= scale;
  last_step *= scale;
  accepted->next_orbital_guess_packed = chart.expand_step(last_step.head(no));
  if (ns) accepted->next_structure_guess =
      structure.coefficient_response_coordinates(last_step.tail(ns));
  record->gradient_dot_step = initial_gradient.dot(displacement);
  record->step_dot_hessian_step = displacement.dot(hessian_displacement);
  record->predicted_reduction =
      -record->gradient_dot_step - 0.5 * record->step_dot_hessian_step;
  record->actual_reduction = energy - accepted->trial.energy;
  record->orbital_step_norm = metric.norm(displacement.head(no));
  record->accepted_trial_radius = record->orbital_step_norm;
  record->next_trust_radius = *trust_radius;
  record->kkt_residual_norm = candidate.residual_norm;
  record->kkt_residual_target = ciah_options.tolerance;
  record->shift = -candidate.eigenvalue;
  // A restored prefix is not the final AH Ritz solution. The recorded AH
  // residual describes that last solve, not a certification of this prefix.
  record->approximate_step = restored_prefix || !candidate.converged;
  return true;
}

}  // namespace xmvb::vb::optimizer_detail
