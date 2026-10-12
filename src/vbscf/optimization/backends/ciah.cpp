#include "vbscf/optimization/backends/ciah.hpp"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <limits>
#include <stdexcept>
#include <utility>

#include <Eigen/Cholesky>
#include "core/metric_solve.hpp"

#include "vbscf/derivatives/hessian/coupled/structure.hpp"
#include "vbscf/derivatives/hessian/exact/operator.hpp"
#include "vbscf/derivatives/hessian/context/accepted_point.hpp"
#include "vbscf/diagnostics/performance/neo_trace.hpp"
#include "vbscf/optimization/driver/options.hpp"
#include "vbscf/optimization/driver/result.hpp"
#include "vbscf/optimization/globalization/line_search.hpp"
#include "vbscf/optimization/neo/ciah.hpp"
#include "vbscf/optimization/preconditioners/coarse.hpp"
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

void CiahMacroStep::align_structure_phase(const Eigen::VectorXd& reference,
    const Eigen::VectorXd& overlap_root) {
  if (reference.size() != overlap_root.size() ||
      next_structure_guess.rows() != reference.size() ||
      !reference.allFinite() || !overlap_root.allFinite())
    throw std::invalid_argument("invalid CIAH restart phase operands");
  if (reference.dot(overlap_root) < 0.0) next_structure_guess *= -1.0;
}

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
  CiahOptions ciah_options;
  double energy_scale = 0.0;
  for (const double state_energy : context->selected_state_energies)
    energy_scale = std::max(energy_scale, std::abs(state_energy));
  ciah_options.operator_relative_accuracy =
      context->structure_solve_accuracy.response_backward_error_tolerance(energy_scale);
  CoarsePreconditioner coarse(dimension,ns ? 2*ciah_options.max_cycle : 0,
      ciah_options.lindep,ciah_options.operator_relative_accuracy);
  Eigen::VectorXd orbital_image(ns ? dimension : 0), structure_image(ns ? dimension : 0);
  record->model_dimension = static_cast<int>(dimension);
  record->source_gradient_norm = gradient.norm();
  record->source_physical_gradient_norm = metric.dual_norm(gradient);

  const auto diagonal_norm = [&](const Eigen::VectorXd& g) {
    return std::hypot(metric.dual_norm(g.head(no)),
        ns ? structure.diagonal_dual_norm(g.tail(ns)) : 0.0);
  };
  const auto dual_norm = [&](const Eigen::VectorXd& g) {
    const bool trace = diagnostics::neo_event_trace_enabled();
    const auto start = trace ? std::chrono::steady_clock::now() :
        std::chrono::steady_clock::time_point{};
    const double orbital_norm = metric.dual_norm(g.head(no));
    const auto structure_start = trace ? std::chrono::steady_clock::now() :
        std::chrono::steady_clock::time_point{};
    const double structure_norm = ns ? structure.dual_norm(g.tail(ns)) : 0.0;
    const auto end = trace ? std::chrono::steady_clock::now() :
        std::chrono::steady_clock::time_point{};
    diagnostics::neo_trace("ciah_coupled_norm", "orbital=%.17g,structure=%.17g",
        orbital_norm, structure_norm);
    if (trace) diagnostics::neo_trace("ciah_norm_cost",
        "orbital_seconds=%.17g,structure_seconds=%.17g",
        std::chrono::duration<double>(structure_start-start).count(),
        std::chrono::duration<double>(end-structure_start).count());
    return std::hypot(orbital_norm, structure_norm);
  };
  // Explicitly bounded diagnostic only; never assemble a production metric.
  // FeCl2 has 419 structure coordinates (about 1.34 MiB per matrix).
  Eigen::LLT<Eigen::MatrixXd> audit_metric;
  const bool audit_structure_metric = ns > 0 &&
      diagnostics::neo_event_trace_enabled() &&
      diagnostics::neo_candidate_audit_enabled() &&
      static_cast<double>(ns) * ns * sizeof(double) <= 16.0 * 1024 * 1024;
  if (audit_structure_metric) {
    Eigen::MatrixXd matrix(ns, ns);
    for (Eigen::Index begin = 0; begin < ns; begin += 32) {
      const Eigen::Index width = std::min<Eigen::Index>(32, ns - begin);
      Eigen::MatrixXd probes = Eigen::MatrixXd::Zero(ns, width);
      probes.middleRows(begin, width).setIdentity();
      matrix.middleCols(begin, width) =
          structure.apply_coupling_coordinate_block(probes).metric_coordinates;
    }
    const double skew = (matrix - matrix.transpose()).norm();
    const double tolerance = 256.0 * std::numeric_limits<double>::epsilon() *
        ns * std::max(1.0, matrix.norm());
    if (skew > tolerance) throw std::runtime_error("CIAH audit metric is not symmetric");
    audit_metric.compute(0.5 * (matrix + matrix.transpose()));
    if (audit_metric.info() != Eigen::Success)
      throw std::runtime_error("CIAH audit metric is not positive definite");
    diagnostics::neo_trace("ciah_structure_metric_audit",
        "dimension=%ld,matrix_bytes=%ld,skew=%.17g",
        static_cast<long>(ns), static_cast<long>(ns * ns * sizeof(double)), skew);
  }
  const auto full_dual_norm = [&](const Eigen::VectorXd& g) {
    const Eigen::VectorXd structure_g = g.tail(ns);
    const Eigen::VectorXd inverse = audit_metric.solve(structure_g);
    const Eigen::VectorXd residual =
        structure.apply_metric_coordinates(inverse) - structure_g;
    diagnostics::neo_trace("ciah_structure_metric_solve_audit",
        "rhs_norm=%.17g,residual=%.17g,solution_norm=%.17g",
        structure_g.norm(), residual.norm(), inverse.norm());
    const double squared = structure_g.dot(inverse);
    if (!(squared >= 0.0) || !std::isfinite(squared))
      throw std::runtime_error("CIAH audit structure dual norm is invalid");
    return std::hypot(metric.dual_norm(g.head(no)), std::sqrt(squared));
  };
  const auto action = [&](const Eigen::VectorXd& x) {
    Eigen::VectorXd image(dimension);
    if (ns == 0) {
      image = hessian.apply_orbital(x);
    } else {
      const auto orbital = hessian.apply_orbital_coupling(x.head(no));
      const auto ci = structure.apply_coupling_coordinates(x.tail(ns));
      const Eigen::VectorXd adjoint = hessian.apply_structure_coupling_adjoint(
          ci.coefficient_response, ci.adjoint_multipliers);
      const Eigen::VectorXd forcing = structure.project_coordinates(
          orbital.scaled_structure_forcing);
      if (orbital.connection.has_value()) {
        const auto& connection = *orbital.connection;
        const Eigen::VectorXd representation = structure.project_coordinates(
            connection.representation_forcing);
        const Eigen::VectorXd integral = structure.project_coordinates(
            connection.integral_forcing);
        const Eigen::VectorXd rhs = structure.project_coordinates(connection.overlap);
        xmvb::core::MetricSolveStats connection_stats;
        const Eigen::VectorXd gamma = xmvb::core::solve_metric(rhs,
            [&](const Eigen::VectorXd& v) { return structure.apply_metric_coordinates(v); },
            [&](const Eigen::VectorXd& v) { return structure.apply_inverse_metric_model(v); },
            &connection_stats);
        const Eigen::VectorXd kinematic = structure.apply_hessian_coordinates(gamma);
        const double representation_norm = structure.dual_norm(representation);
        const double integral_norm = structure.dual_norm(integral);
        const double forcing_norm = structure.dual_norm(forcing);
        diagnostics::neo_trace("ciah_connection",
            "representation=%.17g,integral=%.17g,total=%.17g,cross=%.17g,"
            "reconstruction=%.17g,closure_squared=%.17g,carrier_squared=%.17g,"
            "kinematic_error=%.17g,metric_actions=%d,metric_residual=%.17g",
            representation_norm, integral_norm, forcing_norm,
            0.5 * (forcing_norm * forcing_norm -
                representation_norm * representation_norm - integral_norm * integral_norm),
            (forcing - representation - integral).norm(),
            connection.squared_norm.sum() - rhs.dot(gamma),
            connection.squared_norm.sum(),
            structure.dual_norm(representation - kinematic),
            connection_stats.actions, connection_stats.relative_residual);
      }
      image.head(no) = orbital.orbital_hessian + adjoint;
      image.tail(ns) = forcing + ci.hessian_coordinates;
      orbital_image.head(no) = orbital.orbital_hessian;
      orbital_image.tail(ns) = forcing;
      structure_image.head(no) = adjoint;
      structure_image.tail(ns) = ci.hessian_coordinates;
      if (diagnostics::neo_event_trace_enabled()) {
        const double orbital_energy = x.head(no).dot(orbital.orbital_hessian);
        const double structure_energy = x.tail(ns).dot(ci.hessian_coordinates);
        const double backward_coupling = x.head(no).dot(adjoint);
        const double forward_coupling = x.tail(ns).dot(forcing);
        diagnostics::neo_trace("ciah_block_curvature",
            "orbital=%.17g,structure=%.17g,coupling_forward=%.17g,coupling_backward=%.17g",
            orbital_energy, structure_energy, forward_coupling, backward_coupling);
      }
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
  const auto diagonal_precondition = [&](const Eigen::VectorXd& g, double eigenvalue) {
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
  const auto precondition = [&](const Eigen::VectorXd& g, double eigenvalue) {
    return coarse.apply(g,eigenvalue,[&](const Eigen::VectorXd& v) {
      return diagonal_precondition(v,eigenvalue);
    });
  };
  if (ns) ciah_options.observe_images = [&](const Eigen::VectorXd& x,const Eigen::VectorXd& gx) {
    Eigen::VectorXd v = Eigen::VectorXd::Zero(dimension), g = v;
    v.head(no) = x.head(no);
    g.head(no) = gx.head(no);
    coarse.observe(v,orbital_image,g);
    v.head(no).setZero();
    g.head(no).setZero();
    v.tail(ns) = x.tail(ns);
    g.tail(ns) = gx.tail(ns);
    coarse.observe(v,structure_image,g);
  };
  Eigen::VectorXd initial_gradient = Eigen::VectorXd::Zero(dimension);
  initial_gradient.head(no) = gradient;
  NeoProblem problem(initial_gradient, action, metric_action, dual_norm);
  ciah_options.residual_preconditioner = precondition;
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
  CiahKeyframes keyframes(energy);
  Eigen::VectorXd last_step = displacement;
  double keyframe_norm = dual_norm(initial_gradient);
  double diagonal_keyframe_norm = audit_structure_metric
      ? diagonal_norm(initial_gradient) : 0.0;
  double audit_keyframe_norm = audit_structure_metric
      ? full_dual_norm(initial_gradient) : 0.0;
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
    const double predicted_norm = step_scale == 1.0
        ? candidate.model_gradient_norm : dual_norm(predicted_gradient);
    diagnostics::neo_trace("ciah_micro",
        "basis=%d,micro=%d,omega=%.17g,alpha=%.17g,residual=%.17g,gradient=%.17g,step_scale=%.17g,max_component=%.17g,orbital_residual_euclidean=%.17g,structure_residual_euclidean=%.17g",
        candidate.basis_size, micro_steps, candidate.eigenvalue,
        candidate.reference_weight, candidate.residual_norm, predicted_norm,
        step_scale, largest, candidate.residual.head(no).norm(), candidate.residual.tail(ns).norm());
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
      diagnostics::neo_trace("ciah_micro_stop",
          "reason=%s,micro=%d,budget=%d,basis=%d,keyframe_norm=%.17g,predicted=%.17g,target=%.17g,residual=%.17g,ah_converged=%d",
          micro_stop,micro_steps,micro_limit,candidate.basis_size,keyframe_norm,
          predicted_norm,0.3*micro_gradient_tolerance,candidate.residual_norm,
          static_cast<int>(candidate.converged));
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
      if (diagnostics::neo_event_trace_enabled()) {
        const Eigen::VectorXd correction = exact.gradient - predicted_gradient;
        const double captured = workspace.projected_dual_norm(correction);
        const double secant = correction.dot(keyframe_step);
        diagnostics::neo_trace("ciah_keyframe_secant",
            "mismatch=%.17g,captured=%.17g,displacement_norm=%.17g,secant=%.17g,secant_cosine=%.17g",
            mismatch, captured, length, secant,
            mismatch > 0.0 && length > 0.0 ? secant / (mismatch * length) : 0.0);
      }
      const bool keep_keyframe = mismatch < gradient_trust * predicted_norm ||
          actual_norm < gradient_trust * micro_gradient_tolerance;
      if (audit_structure_metric) {
        const double diagonal_actual = diagonal_norm(exact.gradient);
        const double diagonal_predicted = diagonal_norm(predicted_gradient);
        const double diagonal_mismatch = diagonal_norm(exact.gradient - predicted_gradient);
        const bool diagonal_keep = diagonal_mismatch < gradient_trust * diagonal_predicted ||
            diagonal_actual < gradient_trust * micro_gradient_tolerance;
        const double full_actual = full_dual_norm(exact.gradient);
        const double full_predicted = full_dual_norm(predicted_gradient);
        const double full_mismatch = full_dual_norm(exact.gradient - predicted_gradient);
        const bool full_keep = full_mismatch < gradient_trust * full_predicted ||
            full_actual < gradient_trust * micro_gradient_tolerance;
        diagnostics::neo_trace("ciah_keyframe_metric_audit",
            "diagonal_actual=%.17g,full_actual=%.17g,diagonal_predicted=%.17g,full_predicted=%.17g,diagonal_mismatch=%.17g,full_mismatch=%.17g,diagonal_keep=%d,full_keep=%d,diagonal_trigger=%d,full_trigger=%d",
            diagonal_actual, full_actual, diagonal_predicted, full_predicted,
            diagonal_mismatch, full_mismatch, static_cast<int>(diagonal_keep),
            static_cast<int>(full_keep),
            static_cast<int>(diagonal_predicted < diagonal_keyframe_norm / gradient_trust),
            static_cast<int>(full_predicted < audit_keyframe_norm / gradient_trust));
        if (keep_keyframe) audit_keyframe_norm = full_actual;
      }
      diagnostics::neo_trace("ciah_keyframe",
          "actual=%.17g,predicted=%.17g,mismatch=%.17g,energy=%.17g,accepted=%d",
          actual_norm, predicted_norm, mismatch, exact.energy,
          static_cast<int>(keep_keyframe));
      if (diagnostics::neo_event_trace_enabled()) {
        const double actual_orbital = metric.dual_norm(exact.gradient.head(no));
        const double predicted_orbital = metric.dual_norm(predicted_gradient.head(no));
        const double mismatch_orbital = metric.dual_norm(
            (exact.gradient-predicted_gradient).head(no));
        diagnostics::neo_trace("ciah_keyframe_blocks",
            "actual_orbital=%.17g,actual_structure_squared=%.17g,predicted_orbital=%.17g,predicted_structure_squared=%.17g,mismatch_orbital=%.17g,mismatch_structure_squared=%.17g",
            actual_orbital, std::max(0.0,actual_norm*actual_norm-actual_orbital*actual_orbital),
            predicted_orbital, std::max(0.0,predicted_norm*predicted_norm-predicted_orbital*predicted_orbital),
            mismatch_orbital, std::max(0.0,mismatch*mismatch-mismatch_orbital*mismatch_orbital));
      }
      if (keep_keyframe) {
        keyframes.save(exact.energy,displacement,hessian_displacement,last_step,candidate);
        predicted_gradient = exact.gradient;
        energy_model.refresh(exact.energy, exact.gradient);
        diagnostics::neo_trace("ciah_model_refresh",
            "energy=%.17g,source_energy=%.17g,anchored_prediction=%.17g,origin_prediction=%.17g",
            exact.energy, energy, energy - exact.energy,
            -initial_gradient.dot(displacement) -
                0.5 * displacement.dot(hessian_displacement));
        keyframe_norm = actual_norm;
        if (audit_structure_metric)
          diagonal_keyframe_norm = diagonal_norm(exact.gradient);
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

  double final_reduction = energy - energy_model.energy();
  if (diagnostics::neo_trace_enabled()) {
    // A diagnostic of the actual coupled trial space, not just its A block.
    // Neither inverse callback evaluates an HVP or changes the retained space.
    const auto base = workspace.inverse_spectrum(candidate.eigenvalue,
        [&](const Eigen::VectorXd& v) {
          return diagonal_precondition(v,candidate.eigenvalue);
        });
    const auto balanced = workspace.inverse_spectrum(candidate.eigenvalue,
        [&](const Eigen::VectorXd& v) {
          return precondition(v,candidate.eigenvalue);
        });
    diagnostics::neo_trace("ciah_inverse_spectrum",
        "dimension=%d,positive_curvature=%d,curvature_min=%.17g,curvature_max=%.17g,base_min=%.17g,base_max=%.17g,balanced_min=%.17g,balanced_max=%.17g,base_skew=%.17g,balanced_skew=%.17g",
        base.dimension,static_cast<int>(base.positive_curvature),
        base.curvature_min,base.curvature_max,base.inverse_min,base.inverse_max,
        balanced.inverse_min,balanced.inverse_max,base.inverse_skew,balanced.inverse_skew);
  }
  if (ns && diagnostics::neo_event_trace_enabled()) {
    const auto schur = coarse.schur_spectrum(static_cast<int>(no));
    diagnostics::neo_trace("ciah_schur_spectrum",
        "orbital_dimension=%d,structure_dimension=%d,structure_positive=%d,orbital_positive=%d,"
        "orbital_min=%.17g,orbital_max=%.17g,structure_min=%.17g,structure_max=%.17g,"
        "relaxed_min=%.17g,relaxed_max=%.17g,maximum_fraction_removed=%.17g",
        schur.orbital_dimension,schur.structure_dimension,static_cast<int>(schur.structure_positive),
        static_cast<int>(schur.orbital_positive),schur.orbital_min,schur.orbital_max,
        schur.structure_min,schur.structure_max,schur.relaxed_min,schur.relaxed_max,
        schur.maximum_fraction_removed);
  }
  if (ns && diagnostics::neo_candidate_audit_enabled()) {
    const auto spectrum = coarse.orbital_spectrum(static_cast<int>(no),
        [&](const Eigen::VectorXd& v) { return diagonal_precondition(v,0.0); });
    diagnostics::neo_trace("ciah_orbital_preconditioner_spectrum",
        "dimension=%d,positive_curvature=%d,curvature_min=%.17g,curvature_max=%.17g,preconditioned_min=%.17g,preconditioned_max=%.17g",
        spectrum.dimension, static_cast<int>(spectrum.positive_curvature),
        spectrum.curvature_min, spectrum.curvature_max,
        spectrum.preconditioned_min, spectrum.preconditioned_max);
  }

  accepted->parameters = build_nonredundant_lifted_trial_parameters(
      orbitals, chart, layout, displacement.head(no));
  if ((accepted->parameters - layout.pack(orbitals)).norm() <=
      std::numeric_limits<double>::epsilon() * layout.pack(orbitals).norm() &&
      !keyframes.available()) {
    result->termination_reason = "ciah_zero_macro_step";
    return false;
  }
  // PySCF commits the coupled keyframe sequence and then reoptimizes CI.
  // It does not veto that sequence with the original frozen trust model.
  // Backtrack along evaluated nonlinear keyframes before scaling an origin
  // ray. Only true-energy checkpoints qualify, never obsolete quadratic prefixes.
  double scale = 1.0;
  double linear = initial_gradient.dot(displacement);
  double quadratic = displacement.dot(hessian_displacement);
  double last_ritz_norm = candidate.residual_norm;
  double last_eigenvalue = candidate.eigenvalue;
  bool last_converged = candidate.converged, keyframe_backtracked = false;
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
    if (keyframes.available()) {
      auto point = keyframes.take();
      diagnostics::neo_trace("ciah_keyframe_backtrack",
          "current_energy=%.17g,rejected_energy=%.17g,keyframe_energy=%.17g",
          energy,accepted->trial.energy,point.energy);
      displacement = std::move(point.displacement);
      hessian_displacement = std::move(point.image);
      last_step = std::move(point.last_step);
      last_ritz_norm = point.residual_norm;
      last_eigenvalue = point.eigenvalue;
      last_converged = point.converged;
      final_reduction = energy-point.energy;
      linear = initial_gradient.dot(displacement);
      quadratic = displacement.dot(hessian_displacement);
      scale = 1.0;
      keyframe_backtracked = true;
      accepted->parameters = build_nonredundant_lifted_trial_parameters(
          orbitals,chart,layout,displacement.head(no));
      continue;
    }
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
        last_ritz_norm, scale);
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
  if (ns) {
    // displacement already includes the accepted line-search scale.
    const Eigen::VectorXd reference = coefficients(displacement).col(0);
    const auto& next = *accepted->trial.gradient_result.second_order_context;
    const auto& scf = accepted->trial.gradient_result.scf_result;
    const Eigen::Index roots = next.root_eigenvectors.cols();
    const int root = next.selected_state_indices[0];
    const double raw_overlap = reference.dot(next.selected_state_eigenvectors.col(0));
    Eigen::VectorXd overlap_root;
    if (scf.overlap_eigenvector_matrix.size() ==
        static_cast<std::size_t>(reference.size()*roots)) {
      const Eigen::Map<const Eigen::MatrixXd> mc(
          scf.overlap_eigenvector_matrix.data(),reference.size(),roots);
      overlap_root = mc.col(root);
    } else if (scf.structure_matrices.overlap_matrix.size() ==
        static_cast<std::size_t>(reference.size()*reference.size())) {
      // Dense reference solves already own M; do not assemble it here.
      const Eigen::Map<const Eigen::MatrixXd> m(
          scf.structure_matrices.overlap_matrix.data(),reference.size(),reference.size());
      overlap_root = m*next.selected_state_eigenvectors.col(0);
    } else {
      throw std::logic_error("CIAH restart is missing its relaxed root overlap image");
    }
    accepted->align_structure_phase(reference,overlap_root);
    diagnostics::neo_trace("ciah_seed_phase",
        "raw_overlap=%.17g,metric_overlap=%.17g,seed_norm=%.17g,phase_reversed=%d",
        raw_overlap,reference.dot(overlap_root),accepted->next_structure_guess.norm(),
        static_cast<int>(reference.dot(overlap_root) < 0.0));
  }
  record->gradient_dot_step = initial_gradient.dot(displacement);
  record->step_dot_hessian_step = displacement.dot(hessian_displacement);
  record->predicted_reduction = accepted_prediction;
  record->actual_reduction = energy - accepted->trial.energy;
  record->orbital_step_norm = metric.norm(displacement.head(no));
  record->accepted_trial_radius = record->orbital_step_norm;
  record->next_trust_radius = *trust_radius;
  record->kkt_residual_norm = last_ritz_norm;
  record->kkt_residual_target = ciah_options.tolerance;
  record->shift = -last_eigenvalue;
  // The recorded AH residual belongs to the last microstep, not the complete
  // nonlinear macro displacement.
  record->approximate_step = !last_converged || scale < 1.0 || keyframe_backtracked;
  return true;
}

}  // namespace xmvb::vb::optimizer_detail
