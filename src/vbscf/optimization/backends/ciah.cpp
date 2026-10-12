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
  CiahEnergyModel energy_model(energy, initial_gradient);
  Eigen::VectorXd last_step = displacement;
  double keyframe_norm = dual_norm(initial_gradient);
  // newton_casscf.kernel passes 0.3 * its outer gradient tolerance to
  // update_orb_ci; the microstep/keyframe tests use that inner tolerance.
  const double micro_gradient_tolerance = 0.3 * options.gradient_tolerance;
  int micro_steps = 0, since_keyframe = 0;
  bool audit_next_keyframe_step = false;
  bool audited_keyframe_step = false;
  const char* micro_stop = "unset";
  CiahResult candidate;
  while (true) {
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
    const double step_scale = largest > max_component ? max_component / largest : 1.0;
    if (step_scale < 1.0) {
      step *= step_scale;
      image *= step_scale;
    }
    if (audit_next_keyframe_step) {
      // Check the actual resumed direction at a nonzero keyframe, in the
      // same immutable origin chart used by both the gradient and H images.
      // This opt-in audit neither changes the candidate nor adds any HVP.
      const double length = std::sqrt(std::max(0.0, step.dot(metric_action(step))));
      if (length > 0.0) {
        const Eigen::VectorXd direction = step / length;
        const Eigen::VectorXd origin_image = image / length;
        const double derivative = predicted_gradient.dot(direction);
        for (const double delta : {1.0e-3, 3.0e-4}) {
          const Keyframe plus = evaluate_keyframe(displacement + delta * direction);
          const Keyframe minus = evaluate_keyframe(displacement - delta * direction);
          const double finite = (plus.energy - minus.energy) / (2.0 * delta);
          const Eigen::VectorXd local_image =
              (plus.gradient - minus.gradient) / (2.0 * delta);
          diagnostics::neo_trace("ciah_keyframe_direction_audit",
              "basis=%d,delta=%.17g,analytic=%.17g,finite=%.17g,derivative_error=%.17g,local_curvature=%.17g,origin_curvature=%.17g,curvature_change=%.17g,orbital_derivative=%.17g,structure_derivative=%.17g",
              candidate.basis_size, delta, derivative, finite,
              std::abs(derivative - finite), dual_norm(local_image),
              dual_norm(origin_image), dual_norm(local_image - origin_image),
              predicted_gradient.head(no).dot(direction.head(no)),
              ns ? predicted_gradient.tail(ns).dot(direction.tail(ns)) : 0.0);
        }
      }
      audit_next_keyframe_step = false;
      audited_keyframe_step = true;
    }
    const Eigen::VectorXd previous_gradient = predicted_gradient;
    const Eigen::VectorXd previous_step = last_step;
    displacement += step;
    hessian_displacement += image;
    energy_model.add(step, image);
    predicted_gradient += image;
    last_step = step;
    ++micro_steps;
    ++since_keyframe;
    const double predicted_norm = dual_norm(predicted_gradient);
    diagnostics::neo_trace("ciah_micro",
        "basis=%d,micro=%d,omega=%.17g,alpha=%.17g,residual=%.17g,gradient=%.17g,step_scale=%.17g,max_component=%.17g",
        candidate.basis_size, micro_steps, candidate.eigenvalue,
        candidate.reference_weight, candidate.residual_norm, predicted_norm,
        step_scale, largest);
    if (micro_steps > 3 && predicted_norm > gradient_trust * keyframe_norm) {
      displacement -= step;
      hessian_displacement -= image;
      energy_model.remove(step, image);
      predicted_gradient = previous_gradient;
      last_step = previous_step;
      ++record->rejected_trial_count;
      micro_stop = "gradient_growth";
      break;
    }
    const int micro_limit = std::max(base_micro_cycles,
        base_micro_cycles - static_cast<int>(2.0 * std::log(keyframe_norm + 1.0e-7)));
    if (micro_steps >= micro_limit ||
        predicted_norm < 0.3 * micro_gradient_tolerance || candidate.exhaustion) {
      micro_stop = predicted_norm < 0.3 * micro_gradient_tolerance ? "micro_gradient" :
          candidate.exhaustion ? "subspace_exhaustion" : "micro_budget";
      break;
    }

    // PySCF resets dr after committing a keyframe. The update interval is a
    // function of displacement since that keyframe, not the entire macro.
    const auto& keyframe_step = energy_model.displacement();
    const double length = std::sqrt(std::max(0.0,
        keyframe_step.dot(metric_action(keyframe_step))));
    const double interval = std::max<double>(keyframe_interval,
        keyframe_interval - std::log(length + 1.0e-7));
    if (since_keyframe >= interval || predicted_norm < keyframe_norm / gradient_trust) {
      const Keyframe exact = evaluate_keyframe(displacement);
      const double actual_norm = dual_norm(exact.gradient);
      const double mismatch = dual_norm(exact.gradient - predicted_gradient);
      const bool keep_keyframe = mismatch < gradient_trust * predicted_norm ||
          actual_norm < gradient_trust * micro_gradient_tolerance;
      diagnostics::neo_trace("ciah_keyframe",
          "actual=%.17g,predicted=%.17g,mismatch=%.17g,energy=%.17g,accepted=%d",
          actual_norm, predicted_norm, mismatch, exact.energy,
          static_cast<int>(keep_keyframe));
      if (keep_keyframe) {
        predicted_gradient = exact.gradient;
        energy_model.refresh(exact.energy, exact.gradient);
        diagnostics::neo_trace("ciah_model_refresh",
            "energy=%.17g,source_energy=%.17g,anchored_prediction=%.17g,origin_prediction=%.17g",
            exact.energy, energy, energy - exact.energy,
            -initial_gradient.dot(displacement) -
                0.5 * displacement.dot(hessian_displacement));
        keyframe_norm = actual_norm;
        since_keyframe = 0;
        audit_next_keyframe_step = !audited_keyframe_step &&
            diagnostics::neo_event_trace_enabled() &&
            diagnostics::neo_candidate_audit_enabled();
      } else {
        displacement -= step;
        hessian_displacement -= image;
        energy_model.remove(step, image);
        predicted_gradient = previous_gradient;
        last_step = previous_step;
        ++record->rejected_trial_count;
        micro_stop = "keyframe_mismatch";
        break;
      }
    }
  }

  if (diagnostics::neo_event_trace_enabled()) {
    // These belong to the kept microstep sequence before prefix restoration.
    // A keyframe changes the RHS, not the immutable Hessian images.
    const Eigen::VectorXd frozen_gradient = initial_gradient + hessian_displacement;
    diagnostics::neo_trace("ciah_micro_end",
        "reason=%s,source=%.17g,predicted=%.17g,frozen_model=%.17g,keyframe_correction=%.17g,last_ritz=%.17g",
        micro_stop, dual_norm(initial_gradient), dual_norm(predicted_gradient),
        dual_norm(frozen_gradient), dual_norm(predicted_gradient - frozen_gradient),
        candidate.residual_norm);
  }

  const double final_reduction = energy - energy_model.energy();

  accepted->parameters = build_nonredundant_lifted_trial_parameters(
      orbitals, chart, layout, displacement.head(no));
  if ((accepted->parameters - layout.pack(orbitals)).norm() <=
      std::numeric_limits<double>::epsilon() * layout.pack(orbitals).norm()) {
    result->termination_reason = "ciah_zero_macro_step";
    return false;
  }
  // PySCF commits the coupled keyframe sequence and then reoptimizes CI.
  // It does not veto that sequence with the original frozen trust model.
  // Retain our real-energy descent safeguard, but do not replace the coupled
  // candidate by an old prefix because its obsolete origin model is uphill.
  double scale = 1.0;
  const double linear = initial_gradient.dot(displacement);
  const double quadratic = displacement.dot(hessian_displacement);
  double accepted_prediction = 0.0;
  while (true) {
    accepted->trial = objective->evaluate_trial_with_structure_guess(
        accepted->parameters, coefficients(scale * displacement));
    // Keyframe refreshes both the gradient and the energy anchor. Testing the
    // full microstep sequence against g0*p + p*H0*p/2 discards that update.
    // A shortened origin-chart displacement is a different path: use its
    // origin model, not a scaled prediction for the keyframe sequence.
    const double model_linear = scale == 1.0
        ? -final_reduction - 0.5 * quadratic : scale * linear;
    const double prediction = -model_linear - 0.5 * scale * scale * quadratic;
    const double reduction = energy - accepted->trial.energy;
    const double ratio = prediction > 0.0 ? reduction / prediction
        : -std::numeric_limits<double>::infinity();
    const bool keep_trial = reduction > 0.0;
    diagnostics::neo_trace("ciah_macro_trial",
        "scale=%.17g,linear=%.17g,quadratic=%.17g,current_energy=%.17g,trial_energy=%.17g,predicted=%.17g,actual=%.17g,rho=%.17g,accepted=%d",
        scale, model_linear, scale * scale * quadratic, energy,
        accepted->trial.energy, prediction, reduction, ratio,
        static_cast<int>(keep_trial));
    if (keep_trial) {
      accepted_prediction = prediction;
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
  if (diagnostics::neo_event_trace_enabled()) {
    // r = g0 + H0*p is the unshifted frozen Newton residual of the actual
    // accepted displacement. It reuses accumulated images, adding no HVP.
    // It is NOT the last Ritz residual or a nonlinear stationarity certificate.
    const Eigen::VectorXd frozen_gradient = initial_gradient + hessian_displacement;
    diagnostics::neo_trace("ciah_macro_model",
        "source=%.17g,frozen_model=%.17g,orbital=%.17g,structure_diagonal=%.17g,last_ritz=%.17g,scale=%.17g",
        dual_norm(initial_gradient), dual_norm(frozen_gradient),
        metric.dual_norm(frozen_gradient.head(no)),
        ns ? structure.diagonal_dual_norm(frozen_gradient.tail(ns)) : 0.0,
        candidate.residual_norm, scale);
    if (diagnostics::neo_candidate_audit_enabled()) {
      // One additional fixed-C keyframe, only in explicit nonlinear audit mode.
      // The accepted relaxed gradient is the normal production completion,
      // brought forward and cached so the caller does not repeat its backprop.
      const Keyframe exact = evaluate_keyframe(displacement);
      if (!accepted->trial.gradient_ready) objective->complete_trial(&accepted->trial);
      const Eigen::VectorXd relaxed =
          chart.project_reduced_gradient(accepted->trial.gradient);
      diagnostics::neo_trace("ciah_macro_gradient",
          "source=%.17g,fixed=%.17g,fixed_orbital=%.17g,fixed_structure_diagonal=%.17g,linearization_error=%.17g,relaxed_orbital=%.17g,relaxation_change=%.17g,fixed_energy=%.17g,relaxed_energy=%.17g",
          dual_norm(initial_gradient), dual_norm(exact.gradient),
          metric.dual_norm(exact.gradient.head(no)),
          ns ? structure.diagonal_dual_norm(exact.gradient.tail(ns)) : 0.0,
          dual_norm(exact.gradient - frozen_gradient), metric.dual_norm(relaxed),
          metric.dual_norm(relaxed - exact.gradient.head(no)),
          exact.energy, accepted->trial.energy);
    }
  }
  accepted->next_orbital_guess_packed = chart.expand_step(last_step.head(no));
  if (ns) accepted->next_structure_guess =
      structure.coefficient_response_coordinates(last_step.tail(ns));
  record->gradient_dot_step = initial_gradient.dot(displacement);
  record->step_dot_hessian_step = displacement.dot(hessian_displacement);
  record->predicted_reduction = accepted_prediction;
  record->actual_reduction = energy - accepted->trial.energy;
  record->orbital_step_norm = metric.norm(displacement.head(no));
  record->accepted_trial_radius = record->orbital_step_norm;
  record->next_trust_radius = *trust_radius;
  record->kkt_residual_norm = candidate.residual_norm;
  record->kkt_residual_target = ciah_options.tolerance;
  record->shift = -candidate.eigenvalue;
  // The recorded AH residual belongs to the last microstep, not the complete
  // nonlinear macro displacement.
  record->approximate_step = !candidate.converged || scale < 1.0;
  return true;
}

}  // namespace xmvb::vb::optimizer_detail
