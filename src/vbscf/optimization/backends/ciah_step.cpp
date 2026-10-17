#include <Eigen/SVD>
#include "vbscf/optimization/backends/ciah_step.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <chrono>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

#include "linalg/metric_solve.hpp"

#include "vbscf/derivatives/hessian/coupled/structure.hpp"
#include "vbscf/derivatives/hessian/exact/operator.hpp"
#include "vbscf/derivatives/hessian/context/accepted_point.hpp"
#include "mcscf/optimization/trace.hpp"
#include "vbscf/diagnostics/orbitals/ciah_keyframe.hpp"
#include "vbscf/diagnostics/orbitals/step_geometry.hpp"
#include "core/timing.hpp"
#include "vbscf/optimization/driver/options.hpp"
#include "vbscf/optimization/driver/result.hpp"
#include "vbscf/optimization/globalization/line_search.hpp"
#include "mcscf/optimization/ciah/solver.hpp"
#include "vbscf/optimization/coupled/ci_transport.hpp"
#include "mcscf/optimization/preconditioners/coarse.hpp"
#include "vbscf/optimization/trust_region/retraction.hpp"

namespace xmvb::vb::optimizer_detail {
namespace {

// Initial policies follow PySCF newton_casscf.update_orb_ci. They control
// approximate microsteps, not the VBSCF outer accuracy contract.
constexpr double max_component = 0.03;
constexpr double gradient_trust = 3.0;
constexpr int ah_start_cycle = 3;
constexpr int keyframe_interval = 5;
constexpr int base_micro_cycles = 10;

struct Keyframe {
  Eigen::VectorXd gradient;
  /** Residual-balanced internal energy; macro acceptance remains physical. */
  double energy = 0.0;
  double physical_energy = 0.0;
};

}  // namespace

void CiahMacroStep::align_structure_frame(
    const Eigen::MatrixXd& reference, const Eigen::MatrixXd& overlap_roots) {
  if (reference.rows() != overlap_roots.rows() ||
      reference.cols() != overlap_roots.cols() ||
      next_structure_guess.rows() != reference.rows() ||
      next_structure_guess.cols() != reference.cols())
    throw std::invalid_argument("invalid CIAH restart frame operands");
  const Eigen::MatrixXd cross = reference.transpose() * overlap_roots;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(cross, Eigen::ComputeFullU | Eigen::ComputeFullV);
  next_structure_guess = (next_structure_guess *
      (svd.matrixU() * svd.matrixV().transpose())).eval();
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
    double gradient_norm,
    const Eigen::VectorXd& orbital_guess_packed,
    const Eigen::MatrixXd& structure_guess,
    bool use_structure_response,
    double energy,
    double* trust_radius,
    CiahIterationRecord* record,
    VbScfOptimizerResult* result,
    CiahMacroStep* accepted) {
  xmvb::profiling::ScopedPhase macro_phase(xmvb::profiling::Phase::CiahMacro, true);
  xmvb::profiling::ScopedPhase setup_phase(xmvb::profiling::Phase::CiahSetup);
  const auto context = objective->second_order_context();
  const auto orbitals = objective->input().orbital_preparation_input;
  const Eigen::Index no = chart.reduced_size();
  const Eigen::Index ns = use_structure_response ? structure.tangent_size() : 0;
  const Eigen::Index dimension = no + ns;
  const bool step_audit = diagnostics::orbital_step_audit_enabled();
  std::unique_ptr<CiTransport> transport;
  if (ns) transport = std::make_unique<CiTransport>(
      orbitals, layout, chart, *context, hessian.structure_action(), structure);
  mcscf::CiahOptions ciah_options;
  double energy_scale = 0.0;
  for (const double state_energy : context->selected_state_energies)
    energy_scale = std::max(energy_scale, std::abs(state_energy));
  ciah_options.operator_relative_accuracy =
      context->structure_solve_accuracy.response_backward_error_tolerance(energy_scale);
  mcscf::CoarsePreconditioner coarse(dimension,ns ? 2*ciah_options.max_cycle : 0,
      ciah_options.lindep,ciah_options.operator_relative_accuracy);
  Eigen::VectorXd orbital_image(ns ? dimension : 0), structure_image(ns ? dimension : 0);
  record->model_dimension = static_cast<int>(dimension);
  record->source_gradient_norm = chart.coordinate_dual_norm(gradient);
  record->source_physical_gradient_norm = gradient_norm;

  // Internal CIAH control is in the full coordinate vector, not a projected
  // estimate. Only the cheap local AO inverse is used for coordinate units;
  // there is no complete physical-metric solve for progress/keyframes.
  const auto coordinate_norm = [&](const Eigen::VectorXd& g) {
    return std::hypot(chart.coordinate_dual_norm(g.head(no)), ns ? g.tail(ns).norm() : 0.0);
  };
  const auto action = [&](const Eigen::VectorXd& x, Eigen::VectorXd* metric_image) {
    Eigen::VectorXd image(dimension);
    if (metric_image) {
      metric_image->resize(dimension);
      metric_image->head(no) = metric.apply(x.head(no));
    }
    if (ns == 0) {
      image = hessian.apply_orbital(x);
    } else {
      const bool direct_ci = hessian.structure_action().supports_integral_direction();
      const Eigen::MatrixXd raw_ci=direct_ci ?
          structure.coefficient_response_coordinates(x.tail(ns)):Eigen::MatrixXd{};
      StructureResponseRequest response_request{raw_ci};
      const auto orbital = hessian.apply_orbital_coupling(x.head(no),
          direct_ci ? &response_request:nullptr);
      xmvb::profiling::ScopedPhase structure_phase(xmvb::profiling::Phase::CiahStructure, true);
      const auto structure_start=std::chrono::steady_clock::now();
      StructureActiveIntegralAdjoint active_adjoint;
      const auto ci = direct_ci ? structure.finish_coupling_coordinates(raw_ci,
          std::move(response_request.result),&active_adjoint):
          structure.apply_coupling_coordinates(x.tail(ns));
      const auto coupling_end=std::chrono::steady_clock::now();
      const Eigen::VectorXd adjoint = direct_ci
          ? hessian.apply_structure_coupling_adjoint(active_adjoint)
          : hessian.apply_structure_coupling_adjoint(
              ci.coefficient_response, ci.adjoint_multipliers);
      const auto adjoint_end=std::chrono::steady_clock::now();
      if (metric_image) metric_image->tail(ns) = ci.metric_coordinates;
      const Eigen::VectorXd forcing = structure.project_coordinates(
          orbital.scaled_structure_forcing);
      if (orbital.connection.has_value()) {
        xmvb::profiling::ScopedPhase connection_phase(xmvb::profiling::Phase::CiahConnection, true);
        const auto& connection = *orbital.connection;
        const Eigen::VectorXd representation = structure.project_coordinates(
            connection.representation_forcing);
        const Eigen::VectorXd integral = structure.project_coordinates(
            connection.integral_forcing);
        const Eigen::VectorXd rhs = structure.project_coordinates(connection.overlap);
        xmvb::linalg::MetricSolveStats connection_stats;
        const Eigen::VectorXd gamma = xmvb::linalg::solve_metric(rhs,
            [&](const Eigen::VectorXd& v) { return structure.apply_metric_coordinates(v); },
            [&](const Eigen::VectorXd& v) { return structure.apply_inverse_metric_model(v); },
            &connection_stats);
        const Eigen::VectorXd kinematic = structure.apply_hessian_coordinates(gamma);
        const double representation_norm = structure.dual_norm(representation);
        const double integral_norm = structure.dual_norm(integral);
        const double forcing_norm = structure.dual_norm(forcing);
        mcscf::diagnostics::ciah_trace("ciah_connection",
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
      // The residual-balanced finite model has the physical affine Hessian
      // at this origin, even for a finite CI residual. Its coordinate second
      // jet is cancelled in the keyframe energy AND gradient, not discarded.
      mcscf::diagnostics::ciah_trace("ciah_structure_cost",
          "coupling=%.17g,orbital_adjoint=%.17g,projection=%.17g",
          std::chrono::duration<double>(coupling_end-structure_start).count(),
          std::chrono::duration<double>(adjoint_end-coupling_end).count(),
          std::chrono::duration<double>(std::chrono::steady_clock::now()-adjoint_end).count());
      if (mcscf::diagnostics::ciah_event_trace_enabled()) {
        const double orbital_energy = x.head(no).dot(orbital.orbital_hessian);
        const double structure_energy = x.tail(ns).dot(ci.hessian_coordinates);
        const double backward_coupling = x.head(no).dot(adjoint);
        const double forward_coupling = x.tail(ns).dot(forcing);
        mcscf::diagnostics::ciah_trace("ciah_block_curvature",
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
      direction.head(no) = chart.inverse_coordinate_metric(g.head(no));
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
    xmvb::profiling::ScopedPhase observe_phase(xmvb::profiling::Phase::CiahObserve);
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
  if (ns) initial_gradient.tail(ns) = transport->origin_structure_gradient();
  // QuadraticProblem also serves trust-region solvers that require a physical dual
  // norm. CIAH deliberately does not call this callback.
  mcscf::QuadraticProblem problem(initial_gradient,
      [&](const Eigen::VectorXd& x) { return action(x, nullptr); }, metric_action,
      [&](const Eigen::VectorXd& g) {
        return std::hypot(metric.dual_norm(g.head(no)),
            ns ? structure.dual_norm(g.tail(ns)) : 0.0);
      });
  ciah_options.apply_images = [&](const Eigen::VectorXd& x) {
    Eigen::VectorXd metric_image;
    Eigen::VectorXd hessian_image = action(x, &metric_image);
    return std::make_pair(std::move(hessian_image), std::move(metric_image));
  };
  ciah_options.residual_preconditioner = precondition;
  ciah_options.coordinate_norm = coordinate_norm;
  ciah_options.residual_direction = [&](const Eigen::VectorXd& g) {
    Eigen::VectorXd v = g;
    v.head(no) = chart.inverse_coordinate_metric(g.head(no));
    return v;
  };
  Eigen::VectorXd seed = initial_gradient;
  seed.head(no) = chart.inverse_coordinate_metric(initial_gradient.head(no));
  if (orbital_guess_packed.size() == static_cast<Eigen::Index>(layout.size())) {
    seed.head(no) = chart.project_vector(orbital_guess_packed).reduced_gradient;
  }
  if (ns && structure_guess.rows() == context->selected_state_eigenvectors.rows() &&
      structure_guess.cols() == context->selected_state_eigenvectors.cols()) {
    seed.tail(ns) = structure.project_coordinates(
        std::sqrt(2.0 * context->normalized_state_weights[0]) * structure_guess);
  }
  mcscf::CiahWorkspace workspace(problem, seed, ciah_options);

  const Eigen::MatrixXd c0 = context->selected_state_eigenvectors;
  std::optional<CiTransportPoint> transported;
  Eigen::VectorXd transported_step;
  const auto transport_at = [&](const Eigen::VectorXd& displacement) -> const CiTransportPoint& {
    // One immutable macro model: equal coordinates request the same map.
    // Keyframe and macro-trial consumers share it without another Gram solve.
    if (!transported || transported_step.size()!=displacement.size() ||
        !(transported_step.array()==displacement.array()).all()) {
      transported=transport->retract(displacement.head(no),displacement.tail(ns));
      transported_step=displacement;
    }
    return *transported;
  };
  const auto coefficients = [&](const Eigen::VectorXd& displacement) {
    Eigen::MatrixXd c = c0;
    if (ns) c = transport_at(displacement).coefficients;
    return c;
  };
  const auto evaluate_keyframe = [&](const Eigen::VectorXd& displacement) {
    xmvb::profiling::ScopedPhase keyframe_phase(xmvb::profiling::Phase::CiahKeyframe, true);
    const auto parameters = build_nonredundant_lifted_trial_parameters(
        orbitals, chart, layout, displacement.head(no));
    Keyframe keyframe;
    if (ns) {
      const auto& point = transport_at(displacement);
      const auto exact = objective->evaluate_coupled_keyframe(parameters, point.coefficients);
      auto model = transport->evaluate_model(point, exact);
      keyframe.gradient = std::move(model.gradient);
      keyframe.energy = model.energy;
      keyframe.physical_energy = exact.energy;
    } else {
      const auto exact = objective->evaluate_coupled_keyframe(parameters, c0);
      keyframe.gradient = chart.project_reduced_gradient(exact.orbital_gradient);
      keyframe.energy = exact.energy;
      keyframe.physical_energy = exact.energy;
    }
    ++record->keyframes;
    return keyframe;
  };

  Eigen::VectorXd displacement = Eigen::VectorXd::Zero(dimension);
  Eigen::VectorXd hessian_displacement = displacement;
  Eigen::VectorXd predicted_gradient = initial_gradient;
  mcscf::CiahEnergyModel energy_model(energy, initial_gradient);
  mcscf::CiahKeyframes keyframes(energy);
  Eigen::VectorXd last_step = displacement;
  // Retain the finite Ritz residual rather than assuming exact CI stationarity.
  double keyframe_norm = coordinate_norm(initial_gradient);
  // newton_casscf.kernel passes 0.3 * its outer gradient tolerance to
  // update_orb_ci; the microstep/keyframe tests use that inner tolerance.
  const double micro_gradient_tolerance = 0.3 * options.gradient_tolerance;
  int micro_steps = 0, since_keyframe = 0;
  bool audit_next_keyframe_step = false;
  bool audited_keyframe_step = false;
  bool audited_origin_step = false;
  bool audited_structure_path = false;
  bool energy_model_stopped = false;
  const char* micro_stop = "unset";
  mcscf::CiahResult candidate;
  setup_phase.stop();
  while (true) {
    candidate = workspace.next(predicted_gradient);
    ++record->micro_iterations;
    const double rhs_norm = coordinate_norm(predicted_gradient);
    // A cached-only candidate can be consumed if the CURRENT AH equation is
    // resolved. Otherwise enrich using its current residual before consuming
    // an approximate step. The Ritz residual includes the reference alpha.
    const bool forcing_resolved = candidate.resolved(rhs_norm, 1.0 / gradient_trust);
    const bool ready = candidate.converged || candidate.exhaustion ||
        (candidate.basis_size >= ah_start_cycle &&
         (forcing_resolved || candidate.new_hvp_actions != 0));
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
    if (step_audit) {
      diagnostics::audit_orbital_step(
          orbitals, chart, step.head(no), "micro_origin", micro_steps + 1);
      diagnostics::audit_orbital_step(
          orbitals, chart, step.head(no), "micro_increment", micro_steps + 1,
          displacement.head(no));
    }
    if (!audited_origin_step && displacement.isZero() &&
        mcscf::diagnostics::ciah_trace_enabled() && mcscf::diagnostics::ciah_candidate_audit_enabled()) {
      // Compare the complete derivative in the immutable origin chart. These
      // diagnostic evaluations must not count as optimizer keyframe updates.
      const double length = std::sqrt(std::max(0.0, step.dot(metric_action(step))));
      if (length > 0.0) {
        const Eigen::VectorXd direction = step / length;
        const Eigen::VectorXd exact_image = image / length;
        const double exact_norm = coordinate_norm(exact_image);
        const int previous_keyframes = record->keyframes;
        for (const double delta : {1.0e-3, 3.0e-4}) {
          const Keyframe plus = evaluate_keyframe(delta * direction);
          const Keyframe minus = evaluate_keyframe(-delta * direction);
          const Eigen::VectorXd finite = (plus.gradient - minus.gradient) / (2.0 * delta);
          const Eigen::VectorXd error = finite - exact_image;
          const double derivative = initial_gradient.dot(direction);
          mcscf::diagnostics::ciah_trace("ciah_origin_hvp_audit",
              "delta=%.17g,image_norm=%.17g,finite_norm=%.17g,error_norm=%.17g,orbital_error=%.17g,structure_error=%.17g,derivative=%.17g,finite_derivative=%.17g,orbital_direction_norm=%.17g,structure_direction_norm=%.17g",
              delta, exact_norm, coordinate_norm(finite), coordinate_norm(error),
              error.head(no).norm(),error.tail(ns).norm(),
              derivative, (plus.energy - minus.energy) / (2.0 * delta),
              step.head(no).norm() / length, ns ? step.tail(ns).norm() / length : 0.0);
        }
        record->keyframes = previous_keyframes;
      }
      audited_origin_step = true;
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
          mcscf::diagnostics::ciah_trace("ciah_keyframe_direction_audit",
              "basis=%d,delta=%.17g,analytic=%.17g,finite=%.17g,derivative_error=%.17g,local_curvature=%.17g,origin_curvature=%.17g,curvature_change=%.17g,orbital_derivative=%.17g,structure_derivative=%.17g",
              candidate.basis_size, delta, derivative, finite,
              std::abs(derivative - finite), coordinate_norm(local_image),
              coordinate_norm(origin_image), coordinate_norm(local_image - origin_image),
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
    const double predicted_norm = coordinate_norm(predicted_gradient);
    mcscf::diagnostics::ciah_trace("ciah_micro",
        "basis=%d,micro=%d,omega=%.17g,alpha=%.17g,residual=%.17g,gradient=%.17g,step_scale=%.17g,max_component=%.17g,orbital_residual_euclidean=%.17g,structure_residual_euclidean=%.17g,rhs_norm=%.17g,forcing_target=%.17g,forcing_resolved=%d",
        candidate.basis_size, micro_steps, candidate.eigenvalue,
        candidate.reference_component, candidate.residual_norm, predicted_norm,
        step_scale, largest, candidate.residual.head(no).norm(), candidate.residual.tail(ns).norm(),
        rhs_norm,rhs_norm / gradient_trust,static_cast<int>(forcing_resolved));
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
      mcscf::diagnostics::ciah_trace("ciah_micro_stop",
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
      if (step_audit) diagnostics::audit_orbital_step(
          orbitals, chart, displacement.head(no), "keyframe", micro_steps);
      const Keyframe exact = evaluate_keyframe(displacement);
      if (ns && !audited_structure_path && diagnostics::ciah_keyframe_audit_enabled()) {
        const int saved_keyframes = record->keyframes;
        diagnostics::audit_ciah_keyframe(displacement, hessian_displacement,
            exact.gradient, static_cast<int>(no), [&](const Eigen::VectorXd& point) {
              const auto value = evaluate_keyframe(point);
              return diagnostics::CiahKeyframeProbe{value.gradient, value.energy};
            });
        record->keyframes = saved_keyframes;
        audited_structure_path = true;
      }
      const double actual_norm = coordinate_norm(exact.gradient);
      const Eigen::VectorXd correction = exact.gradient - predicted_gradient;
      const double mismatch_norm = coordinate_norm(correction);
      const bool small_gradient = actual_norm < gradient_trust * micro_gradient_tolerance;
      const bool mismatch_small = mismatch_norm < gradient_trust * predicted_norm;
      if (mcscf::diagnostics::ciah_event_trace_enabled()) {
        const double captured = workspace.projected_dual_norm(correction);
        const double secant = correction.dot(keyframe_step);
        mcscf::diagnostics::ciah_trace("ciah_keyframe_secant",
            "mismatch_norm=%.17g,captured=%.17g,displacement_norm=%.17g,secant=%.17g,mismatch_small=%d",
            mismatch_norm, captured, length, secant, static_cast<int>(mismatch_small));
      }
      const bool keep_keyframe = small_gradient || mismatch_small;
      // Model agreement controls further microsteps, not true-energy descent.
      // Keep an evaluated lower-energy candidate even when its model is stale;
      // normal macro acceptance still reoptimizes CI and checks real energy.
      keyframes.save(exact.energy,displacement,hessian_displacement,last_step,candidate);
      mcscf::diagnostics::ciah_trace("ciah_keyframe",
          "actual=%.17g,predicted=%.17g,mismatch=%.17g,mismatch_small=%d,energy=%.17g,physical_energy=%.17g,accepted=%d",
          actual_norm, predicted_norm, mismatch_norm, static_cast<int>(mismatch_small),
          exact.energy, exact.physical_energy, static_cast<int>(keep_keyframe));
      if (mcscf::diagnostics::ciah_event_trace_enabled()) {
        const double actual_orbital = chart.coordinate_dual_norm(exact.gradient.head(no));
        const double predicted_orbital = chart.coordinate_dual_norm(predicted_gradient.head(no));
        const double mismatch_orbital = chart.coordinate_dual_norm(correction.head(no));
        mcscf::diagnostics::ciah_trace("ciah_keyframe_blocks",
            "actual_orbital=%.17g,actual_structure_squared=%.17g,predicted_orbital=%.17g,predicted_structure_squared=%.17g,mismatch_orbital=%.17g,mismatch_structure_squared=%.17g",
            actual_orbital, std::max(0.0,actual_norm*actual_norm-actual_orbital*actual_orbital),
            predicted_orbital, std::max(0.0,predicted_norm*predicted_norm-predicted_orbital*predicted_orbital),
            mismatch_orbital, std::max(0.0,mismatch_norm*mismatch_norm-mismatch_orbital*mismatch_orbital));
      }
      if (keep_keyframe) {
        // Gradient agreement alone can miss an energy model that predicts
        // descent while its exact segment climbs. Do not extrapolate that
        // frozen curvature further. Keep this actual keyframe displacement
        // for the normal relaxed-CI macro trial; it may still lower energy.
        if (energy_model.descent_mismatch(exact.energy)) {
          energy_model_stopped = true;
          micro_stop = "keyframe_energy_model";
          mcscf::diagnostics::ciah_trace("ciah_keyframe_energy_model",
              "anchor=%.17g,predicted=%.17g,actual=%.17g,stop=1",
              energy_model.anchor_energy(),energy_model.energy(),exact.energy);
          break;
        }
        xmvb::profiling::ScopedPhase refresh_phase(xmvb::profiling::Phase::CiahRefresh);
        predicted_gradient = exact.gradient;
        energy_model.refresh(exact.energy, exact.gradient);
        mcscf::diagnostics::ciah_trace("ciah_model_refresh",
            "energy=%.17g,source_energy=%.17g,anchored_prediction=%.17g,origin_prediction=%.17g",
            exact.energy, energy, energy - exact.energy,
            -initial_gradient.dot(displacement) -
                0.5 * displacement.dot(hessian_displacement));
        keyframe_norm = actual_norm;
        since_keyframe = 0;
        audit_next_keyframe_step = !audited_keyframe_step &&
            mcscf::diagnostics::ciah_event_trace_enabled() &&
            mcscf::diagnostics::ciah_candidate_audit_enabled();
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
    // The microstep and possibly a keyframe changed the RHS. Reproject that
    // RHS in the retained space before spending an HVP on enrichment; a
    // residual from the previously consumed equation is no longer current.
    workspace.reset_rhs();
  }

  if (mcscf::diagnostics::ciah_event_trace_enabled()) {
    // These belong to the kept microstep sequence before prefix restoration.
    // A keyframe changes the RHS, not the immutable Hessian images.
    const Eigen::VectorXd frozen_gradient = initial_gradient + hessian_displacement;
    mcscf::diagnostics::ciah_trace("ciah_micro_end",
        "reason=%s,source_coordinate_norm=%.17g,predicted_coordinate_norm=%.17g,frozen_model_coordinate_norm=%.17g,keyframe_correction_norm=%.17g,last_ritz=%.17g",
        micro_stop, coordinate_norm(initial_gradient), coordinate_norm(predicted_gradient),
        coordinate_norm(frozen_gradient), coordinate_norm(predicted_gradient - frozen_gradient),
        candidate.residual_norm);
  }

  double final_reduction = energy - energy_model.energy();
  if (mcscf::diagnostics::ciah_trace_enabled()) {
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
    mcscf::diagnostics::ciah_trace("ciah_inverse_spectrum",
        "dimension=%d,positive_curvature=%d,curvature_min=%.17g,curvature_max=%.17g,base_min=%.17g,base_max=%.17g,balanced_min=%.17g,balanced_max=%.17g,base_skew=%.17g,balanced_skew=%.17g",
        base.dimension,static_cast<int>(base.positive_curvature),
        base.curvature_min,base.curvature_max,base.inverse_min,base.inverse_max,
        balanced.inverse_min,balanced.inverse_max,base.inverse_skew,balanced.inverse_skew);
  }
  if (ns && mcscf::diagnostics::ciah_event_trace_enabled()) {
    const auto schur = coarse.schur_spectrum(static_cast<int>(no));
    mcscf::diagnostics::ciah_trace("ciah_schur_spectrum",
        "orbital_dimension=%d,structure_dimension=%d,structure_positive=%d,orbital_positive=%d,"
        "orbital_min=%.17g,orbital_max=%.17g,structure_min=%.17g,structure_max=%.17g,"
        "relaxed_min=%.17g,relaxed_max=%.17g,maximum_fraction_removed=%.17g",
        schur.orbital_dimension,schur.structure_dimension,static_cast<int>(schur.structure_positive),
        static_cast<int>(schur.orbital_positive),schur.orbital_min,schur.orbital_max,
        schur.structure_min,schur.structure_max,schur.relaxed_min,schur.relaxed_max,
        schur.maximum_fraction_removed);
  }
  if (ns && mcscf::diagnostics::ciah_candidate_audit_enabled()) {
    const auto spectrum = coarse.orbital_spectrum(static_cast<int>(no),
        [&](const Eigen::VectorXd& v) { return diagonal_precondition(v,0.0); });
    mcscf::diagnostics::ciah_trace("ciah_orbital_preconditioner_spectrum",
        "dimension=%d,positive_curvature=%d,curvature_min=%.17g,curvature_max=%.17g,preconditioned_min=%.17g,preconditioned_max=%.17g",
        spectrum.dimension, static_cast<int>(spectrum.positive_curvature),
        spectrum.curvature_min, spectrum.curvature_max,
        spectrum.preconditioned_min, spectrum.preconditioned_max);
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
  bool origin_restarted = false;
  // Keyframe corrections can turn the accumulated path into an uphill ray
  // from the macro origin. Shrinking that ray cannot restore local descent.
  // Re-solve the ORIGIN augmented model in the already retained AH space;
  // its lowest root supplies descent or negative curvature without a new HVP.
  const auto restart_origin = [&]() {
    if (origin_restarted) return false;
    const auto ray = workspace.origin_ray(initial_gradient);
    if (!ray) return false;
    double largest = chart.maximum_rotation_component(ray->step.head(no));
    if (ns) largest = std::max(largest,
        structure.maximum_coefficient_component_coordinates(ray->step.tail(ns)));
    if (!(largest > 0.0)) return false;
    // Clip BEFORE dividing by a possibly tiny AH reference component. A
    // Hessian-only negative-curvature root has alpha=0 and uses the same cap.
    const double factor = largest > max_component * ray->reference_component
        ? max_component / largest : 1.0 / ray->reference_component;
    displacement = factor * ray->step;
    hessian_displacement = factor * ray->hessian_step;
    last_step = displacement;
    linear = initial_gradient.dot(displacement);
    quadratic = displacement.dot(hessian_displacement);
    final_reduction = -linear - 0.5 * quadratic;
    last_ritz_norm = ray->residual_norm;
    last_eigenvalue = ray->eigenvalue;
    last_converged = false;
    scale = 1.0;
    origin_restarted = true;
    accepted->parameters = build_nonredundant_lifted_trial_parameters(
        orbitals, chart, layout, displacement.head(no));
    mcscf::diagnostics::ciah_trace("ciah_origin_restart",
        "basis=%d,linear=%.17g,quadratic=%.17g,alpha=%.17g,eigenvalue=%.17g,residual=%.17g,hvp_added=0",
        workspace.basis_size(),linear,quadratic,ray->reference_component,
        ray->eigenvalue,ray->residual_norm);
    return true;
  };
  accepted->parameters = build_nonredundant_lifted_trial_parameters(
      orbitals, chart, layout, displacement.head(no));
  if ((accepted->parameters - layout.pack(orbitals)).norm() <=
      std::numeric_limits<double>::epsilon() * layout.pack(orbitals).norm() &&
      !keyframes.available() && !restart_origin()) {
    result->termination_reason = "ciah_zero_macro_step";
    return false;
  }
  double accepted_prediction = 0.0;
  while (true) {
    xmvb::profiling::ScopedPhase trial_phase(xmvb::profiling::Phase::CiahTrial, true);
    if (step_audit) diagnostics::audit_orbital_step(
        orbitals, chart, scale * displacement.head(no), "macro_trial", micro_steps);
    accepted->trial = objective->evaluate_trial_energy_with_structure_guess(
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
    mcscf::diagnostics::ciah_trace("ciah_macro_trial",
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
      mcscf::diagnostics::ciah_trace("ciah_keyframe_backtrack",
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
    if ((energy_model_stopped || linear >= 0.0) && restart_origin()) continue;
    scale *= 0.5;
    if (scale * metric.norm(displacement.head(no)) <= options.minimum_step_size) {
      if (restart_origin()) continue;
      // Failure evidence from existing quantities only: no gradient, HVP or
      // metric inverse is added, and the acceptance policy is unchanged.
      std::fprintf(stderr,
          "CIAH_FAILURE,iteration=%d,reason=ciah_macro_energy_rejected,"
          "trials=%d,origin_slope=%.17g,origin_curvature=%.17g,"
          "source_energy=%.17g,last_trial_energy=%.17g,next_scale=%.17g,"
          "keyframe_backtracked=%d,micro_stop=%s,hvp=%d\n",
          record->accepted_iteration_index,record->rejected_trial_count,
          linear,quadratic,energy,accepted->trial.energy,scale,
          static_cast<int>(keyframe_backtracked),micro_stop,record->orbital_hvp_actions);
      result->termination_reason = "ciah_macro_energy_rejected";
      return false;
    }
    accepted->parameters = build_nonredundant_lifted_trial_parameters(
        orbitals, chart, layout, scale * displacement.head(no));
  }
  // Only the accepted point needs gradient backpropagation and its context.
  // Complete the retained forward evaluation, rather than rebuilding it.
  xmvb::profiling::ScopedPhase completion_phase(xmvb::profiling::Phase::CiahCompletion, true);
  objective->complete_trial(&accepted->trial);
  completion_phase.stop();
  displacement *= scale;
  hessian_displacement *= scale;
  last_step *= scale;
  if (step_audit) diagnostics::audit_orbital_step(
      orbitals, chart, displacement.head(no), "macro_accepted", micro_steps);
  if (mcscf::diagnostics::ciah_event_trace_enabled()) {
    // r = g0 + H0*p is the unshifted frozen Newton residual of the actual
    // accepted displacement. It reuses accumulated images, adding no HVP.
    // It is NOT the last Ritz residual or a nonlinear stationarity certificate.
    const Eigen::VectorXd frozen_gradient = initial_gradient + hessian_displacement;
    mcscf::diagnostics::ciah_trace("ciah_macro_model",
        "source_coordinate_norm=%.17g,frozen_model_coordinate_norm=%.17g,orbital=%.17g,structure_coordinate_norm=%.17g,last_ritz=%.17g,scale=%.17g",
        coordinate_norm(initial_gradient), coordinate_norm(frozen_gradient),
        frozen_gradient.head(no).norm(),
        frozen_gradient.tail(ns).norm(),
        last_ritz_norm, scale);
    if (mcscf::diagnostics::ciah_candidate_audit_enabled()) {
      // One additional fixed-C keyframe, only in explicit nonlinear audit mode.
      // The accepted relaxed gradient is the normal production completion,
      // brought forward and cached so the caller does not repeat its backprop.
      const Keyframe exact = evaluate_keyframe(displacement);
      if (!accepted->trial.gradient_ready) objective->complete_trial(&accepted->trial);
      const Eigen::VectorXd relaxed =
          chart.project_reduced_gradient(accepted->trial.gradient);
      mcscf::diagnostics::ciah_trace("ciah_macro_gradient",
          "source=%.17g,fixed=%.17g,fixed_orbital=%.17g,fixed_structure_coordinate_norm=%.17g,linearization_error=%.17g,relaxed_orbital=%.17g,relaxation_change=%.17g,fixed_energy=%.17g,relaxed_energy=%.17g",
          coordinate_norm(initial_gradient), coordinate_norm(exact.gradient),
          exact.gradient.head(no).norm(),
          exact.gradient.tail(ns).norm(),
          coordinate_norm(exact.gradient - frozen_gradient), relaxed.norm(),
          (relaxed - exact.gradient.head(no)).norm(),
          exact.energy, accepted->trial.energy);
    }
  }
  accepted->next_orbital_guess_packed = chart.expand_step(last_step.head(no));
  if (ns) accepted->next_structure_guess =
      structure.coefficient_response_coordinates(last_step.tail(ns));
  if (ns) {
    const Eigen::MatrixXd reference = coefficients(displacement);
    const auto& next = *accepted->trial.gradient_result.second_order_context;
    accepted->align_structure_frame(reference,
        next.selected_state_structure_images.value().overlap);
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
  // This is the last microstep Ritz residual, or the unscaled origin-root
  // residual after restart. Neither certifies the accepted nonlinear step.
  record->approximate_step = !last_converged || scale < 1.0 ||
      keyframe_backtracked || origin_restarted || energy_model_stopped;
  return true;
}

}  // namespace xmvb::vb::optimizer_detail
