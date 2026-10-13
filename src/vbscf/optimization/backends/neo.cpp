#include "vbscf/optimization/backends/neo.hpp"
#include "vbscf/optimization/backends/ciah.hpp"
#include "vbscf/diagnostics/performance/neo_trace.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <functional>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <Eigen/Eigenvalues>

#include "vbscf/derivatives/hessian/exact/operator.hpp"
#include "vbscf/derivatives/hessian/context/accepted_point.hpp"
#include "vbscf/derivatives/hessian/coupled/structure.hpp"
#include "vbscf/optimization/driver/checks.hpp"
#include "vbscf/optimization/driver/options.hpp"
#include "vbscf/optimization/driver/result.hpp"
#include "vbscf/optimization/driver/session.hpp"
#include "vbscf/optimization/globalization/line_search.hpp"
#include "vbscf/optimization/neo/globalization.hpp"
#include "vbscf/optimization/neo/response_solver.hpp"
#include "vbscf/optimization/objective/function.hpp"
#include "vbscf/optimization/preconditioners/shifted_metric.hpp"
#include "vbscf/optimization/trust_region/retraction.hpp"
#include "vbscf/orbitals/charts/layout.hpp"
#include "vbscf/structures/assembly/selected_coefficients.hpp"

namespace xmvb::vb::optimizer_detail {
namespace {

OrbitalPreconditioner resolve_neo_preconditioner(
    const VbScfObjective& objective,
    OrbitalPreconditioner requested) {
  if (requested != OrbitalPreconditioner::Automatic) return requested;
  return objective.second_order_context()->selected_state_indices.size() == 1
      ? OrbitalPreconditioner::HessianDiagonal
      : OrbitalPreconditioner::OneElectron;
}

double neo_forcing_term(
    double gradient_norm,
    Eigen::Index dimension) {
  // eta = ||g_k||/(1+||g_k||) is below one globally and is O(||g_k||)
  // locally, as required for an inexact Newton method to retain quadratic
  // convergence. Unlike normalization by the initial gradient, it remains
  // meaningful when an initial guess happens to lie near a saddle point.
  // The lower bound only prevents asking for accuracy below accumulated
  // floating-point roundoff.
  const double forcing = gradient_norm / (1.0 + gradient_norm);
  const double roundoff = std::numeric_limits<double>::epsilon() *
      static_cast<double>(std::max<Eigen::Index>(1, dimension));
  return std::max(roundoff, forcing);
}

bool admit_structure_response(
    const VbScfObjective& objective,
    const StructureTangentOperator& structure_hessian,
    const StructureAction& structure_action,
    Eigen::Index orbital_dimension) {
  if (structure_hessian.tangent_size() == 0) return false;

  const auto& context = *objective.second_order_context();
  const auto& orbitals = objective.input().orbital_preparation_input;
  const long double n_ao = orbitals.n_basis_functions;
  const long double n_active = orbitals.n_active_orbitals;
  const long double n_ao_pairs = n_ao * (n_ao + 1.0L) / 2.0L;
  const long double n_active_pairs =
      n_active * (n_active + 1.0L) / 2.0L;
  const auto& graph = objective.input().ao_integral_input.pair_graph;
  // Exact HVP applies the AO-pair operator to active-pair columns. Its work
  // scales with graph edges, not just the number of AO-pair rows.
  const long double ao_pair_actions = graph.empty()
      ? n_ao_pairs
      : graph.has_row_graph()
          ? static_cast<long double>(graph.values.size())
          : 2.0L * static_cast<long double>(graph.integral_count());
  const long double core_pair_work = ao_pair_actions * n_active_pairs;
  // The orthogonal direct-CI branch uses pair-channel GEMM on spin products,
  // not the U_alpha^2 U_beta + U_alpha U_beta^2 cofactor sandwiches.
  const auto& selected = context.selected_state_matrices;
  const long double response_work =
      structure_action.supports_integral_direction()
          ? static_cast<long double>(selected.n_unique_alpha) *
              selected.n_unique_beta * n_active_pairs * n_active_pairs *
              selected.states.size()
          : estimate_selected_state_contraction_work(selected);

  // A coupled orbital basis vector evaluates Bp in addition to the core HVP.
  // Retain the complete coupled model when either its horizontal response
  // space does not exceed the orbital tangent space or its leading
  // unique-string contraction is no more expensive than the core AO-pair
  // contraction.  Otherwise use orbital-only NEO and let the exact trial
  // energy and gradient globalize it.  Both tests compare mathematical work
  // dimensions; neither contains a molecule, active-space, or timing cutoff.
  return structure_hessian.tangent_size() <= orbital_dimension ||
      response_work <= core_pair_work;
}

struct AcceptedNeoKeyframe {
  VbScfObjective::TrialEvaluation trial;
  Eigen::VectorXd parameters;
  Eigen::VectorXd next_orbital_guess_packed;
  Eigen::MatrixXd next_structure_guess;
};

/** @brief Measures a completed trial gradient in its own physical chart. */
double trial_physical_gradient_norm(
    const VbScfObjective::TrialEvaluation& trial,
    const SparseParameterLayout& parameter_view,
    bool complete_active_space) {
  const auto& orbitals = trial.orbital_preparation_input;
  const auto& preparation = trial.gradient_result.orbital_preparation_result;
  const int n_inactive = (orbitals.n_total_electrons - orbitals.n_active_electrons) / 2;
  const OrbitalChart chart(
      orbitals, parameter_view,
      preparation.auxiliary_orbital_matrix.leftCols(n_inactive + orbitals.n_active_orbitals),
      preparation.physical_orbital_frame.normalized_orbital_matrix,
      nullptr, false, complete_active_space, nullptr, OrbitalChart::Coordinates::Unwhitened);
  const NonredundantRetractionMetric metric(chart, parameter_view, orbitals);
  return metric.dual_norm(chart.project_gradient(trial.gradient).reduced_gradient);
}

bool build_accepted_neo_keyframe(
    VbScfObjective* objective,
    const SparseParameterLayout& parameter_view,
    const VbScfOptimizerOptions& options,
    const OrbitalChart& chart,
    const Eigen::VectorXd& reduced_gradient,
    const Eigen::VectorXd& current_parameters,
    const Eigen::VectorXd& orbital_guess_packed,
    const Eigen::MatrixXd& structure_guess,
    double energy,
    double* trust_radius,
    NeoIterationRecord* record,
    VbScfOptimizerResult* result,
    AcceptedNeoKeyframe* accepted) {
  const OrbitalPreparationInput accepted_orbitals =
      objective->input().orbital_preparation_input;
  NonredundantRetractionMetric orbital_metric(
      chart, parameter_view, accepted_orbitals);
  record->source_gradient_norm = chart.coordinate_dual_norm(reduced_gradient);
  const double physical_gradient_norm = orbital_metric.dual_norm(reduced_gradient);
  record->source_physical_gradient_norm = physical_gradient_norm;
  ExactHvpOperator orbital_hessian(
      objective->second_order_context(),
      &objective->input(),
      parameter_view,
      &chart);
  if (!orbital_hessian.supports_analytic_core_model()) {
    result->termination_reason = "neo_analytic_hessian_unavailable";
    return false;
  }
  StructureTangentOperator structure_hessian(
      objective->second_order_context(), orbital_hessian.structure_action());
  const bool use_structure_response =
      admit_structure_response(
          *objective, structure_hessian, orbital_hessian.structure_action(),
          chart.reduced_size());
  if (objective->second_order_context()->selected_state_indices.size() == 1) {
    diagnostics::neo_trace_context =
        {record->accepted_iteration_index, 0, "ciah"};
    CiahMacroStep macro;
    if (!build_ciah_macro_step(objective, parameter_view, options, chart,
            orbital_metric, orbital_hessian, structure_hessian,
            reduced_gradient, physical_gradient_norm, orbital_guess_packed, structure_guess, use_structure_response,
            energy, trust_radius, record, result, &macro)) return false;
    accepted->trial = std::move(macro.trial);
    accepted->parameters = std::move(macro.parameters);
    accepted->next_orbital_guess_packed = std::move(macro.next_orbital_guess_packed);
    accepted->next_structure_guess = std::move(macro.next_structure_guess);
    return true;
  }
  record->model_dimension = std::max(
      record->model_dimension,
      static_cast<int>(
          chart.reduced_size() +
          (use_structure_response ? structure_hessian.tangent_size() : 0)));
  Eigen::VectorXd orbital_guess;
  if (orbital_guess_packed.size() ==
      static_cast<Eigen::Index>(parameter_view.size())) {
    orbital_guess =
        chart.project_vector(orbital_guess_packed).reduced_gradient;
  }
  double selected_energy_scale = 0.0;
  for (const double selected_energy :
       objective->second_order_context()->selected_state_energies) {
    selected_energy_scale =
        std::max(selected_energy_scale, std::abs(selected_energy));
  }
  const double operator_relative_accuracy =
      objective->second_order_context()
          ->structure_solve_accuracy
          .response_backward_error_tolerance(selected_energy_scale);
  // Evaluate one candidate at a time, retaining no direction or trial cache.
  // Counters are per radius attempt, matching ResponseNeoResult's HVP count.
  std::function<void(const ResponseNeoResult&)> observe_candidate;
  if (diagnostics::neo_candidate_audit_enabled()) {
    observe_candidate = [&, last_attempt = -1, last_hvp = -1](
                            const ResponseNeoResult& candidate) mutable {
      const int attempt = diagnostics::neo_trace_context.attempt;
      if (attempt != last_attempt) {
        last_attempt = attempt;
        last_hvp = -1;
      }
      const int hvp = candidate.orbital_actions;
      if (hvp == last_hvp ||
          (hvp != 3 && hvp != 5 && hvp != 10 && hvp != 20 &&
           hvp != 50 && hvp != 100)) return;
      last_hvp = hvp;
      const auto start = std::chrono::steady_clock::now();
      const auto parameters = build_nonredundant_lifted_trial_parameters(
          accepted_orbitals, chart, parameter_view, candidate.step.orbital);
      // Trial evaluation restores the accepted input and is never committed.
      auto trial = objective->evaluate_trial(parameters, true);
      const double actual_gradient = trial_physical_gradient_norm(
          trial, parameter_view, objective->input().complete_active_space);
      const double predicted_gradient = orbital_metric.dual_norm(
          reduced_gradient + candidate.hessian_step.orbital);
      const double gradient_dot_step = reduced_gradient.dot(candidate.step.orbital);
      const double curvature = candidate.step.orbital.dot(candidate.hessian_step.orbital) +
          candidate.step.structure.dot(candidate.hessian_step.structure);
      const auto globalization = globalize_neo_trial(
          energy, trial.energy, gradient_dot_step, curvature,
          candidate.boundary, *trust_radius, options.minimum_step_size);
      const double seconds = std::chrono::duration<double>(
          std::chrono::steady_clock::now() - start).count();
      diagnostics::neo_trace("quality",
          "hvp=%d,micro=%d,radius=%.17g,boundary=%d,shift=%.17g,residual=%.17g,target=%.17g,predicted=%.17g,actual=%.17g,rho=%.17g,energy_acceptable=%d,source_gradient=%.17g,actual_gradient=%.17g,predicted_gradient=%.17g,audit_seconds=%.17g",
          hvp, candidate.iterations, *trust_radius, static_cast<int>(candidate.boundary),
          candidate.shift, candidate.residual_norm, candidate.residual_target,
          globalization.predicted_reduction, globalization.actual_reduction,
          globalization.rho, static_cast<int>(globalization.accepted),
          physical_gradient_norm, actual_gradient, predicted_gradient, seconds);
    };
  }
  if (diagnostics::neo_trace_enabled()) {
    const int n_inactive =
        (accepted_orbitals.n_total_electrons - accepted_orbitals.n_active_electrons) / 2;
    observe_candidate = [&, n_inactive, audit = std::move(observe_candidate)](
                            const ResponseNeoResult& candidate) {
      if (audit) audit(candidate);
      const int count = candidate.orbital_actions;
      if (count <= 0 || ((count & (count - 1)) != 0 &&
                        candidate.residual_norm > candidate.residual_target)) return;
      const Eigen::VectorXd& residual = candidate.kkt_residual.orbital;
      const Eigen::VectorXd inverse_residual =
          chart.apply_inverse_reduced_block_preconditioner(residual);
      const Eigen::VectorXd weights = chart.orbital_pairings(inverse_residual, residual);
      const double total = weights.sum();
      for (Eigen::Index orbital = 0; orbital < weights.size(); ++orbital) {
        diagnostics::neo_trace("orbital_residual",
            "hvp=%d,orbital=%ld,inactive=%d,block_pairing=%.17g,fraction=%.17g,total=%.17g",
            count, static_cast<long>(orbital + 1), static_cast<int>(orbital < n_inactive),
            weights[orbital], total > 0.0 ? weights[orbital] / total : 0.0, total);
      }
    };
  }
  std::function<ResponseNeoResult(const NeoOptions&, bool)> solve_step;
  if (use_structure_response) {
    ResponseNeoModelObserver observe_model;
    if (diagnostics::neo_trace_enabled()) {
      observe_model = [&chart, last_dimension = Eigen::Index{0}](
                          const Eigen::MatrixXd& basis,
                          const Eigen::MatrixXd& core,
                          const Eigen::MatrixXd& model) mutable {
        // Optional diagnostic export only. Do not retain or duplicate the
        // large-space basis in the production solver.
        if (const char* prefix = std::getenv("XMVB_NEO_MODEL_AUDIT")) {
          const auto write = [prefix](const char* suffix, const Eigen::MatrixXd& matrix) {
            std::ofstream stream(std::string(prefix) + suffix, std::ios::binary);
            const std::int64_t shape[] = {matrix.rows(), matrix.cols()};
            stream.write(reinterpret_cast<const char*>(shape), sizeof(shape));
            stream.write(reinterpret_cast<const char*>(matrix.data()),
                         matrix.size() * sizeof(double));
            if (!stream) throw std::runtime_error("failed to export NEO model audit");
          };
          write("-basis.bin", basis);
          write("-model.bin", model);
        }
        const Eigen::Index dimension = basis.cols();
        // Sample doubling dimensions only; structure enrichments at the same
        // dimension must not trigger repeated diagnostic matrix products.
        if (dimension == last_dimension || (dimension & (dimension - 1)) != 0) return;
        last_dimension = dimension;
        Eigen::MatrixXd images(basis.rows(), dimension);
        for (Eigen::Index column = 0; column < dimension; ++column)
          images.col(column) = chart.apply_reduced_curvature(basis.col(column));
        const Eigen::MatrixXd projected = basis.transpose() * images;
        const Eigen::MatrixXd symmetric = 0.5 * (projected + projected.transpose());
        diagnostics::neo_trace("model_difference",
            "dimension=%ld,block=%.17g,fixed=%.17g,response=%.17g,fixed_defect=%.17g",
            static_cast<long>(dimension), symmetric.stableNorm(), core.stableNorm(),
            (model - core).stableNorm(), (core - symmetric).stableNorm());
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(
            model, symmetric);
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> core_spectrum(
            core, symmetric, Eigen::EigenvaluesOnly);
        if (spectrum.info() != Eigen::Success || core_spectrum.info() != Eigen::Success)
          throw std::runtime_error("failed to audit NEO block-preconditioned spectrum");
        diagnostics::neo_trace("preconditioned_spectrum",
            "dimension=%ld,min=%.17g,max=%.17g,core_min=%.17g,core_max=%.17g",
            static_cast<long>(dimension), spectrum.eigenvalues()[0],
            spectrum.eigenvalues()[dimension - 1], core_spectrum.eigenvalues()[0],
            core_spectrum.eigenvalues()[dimension - 1]);
        const Eigen::VectorXd soft = basis * spectrum.eigenvectors().col(0);
        const Eigen::VectorXd soft_image = images * spectrum.eigenvectors().col(0);
        const Eigen::VectorXd weights = chart.orbital_pairings(soft, soft_image);
        for (Eigen::Index orbital = 0; orbital < weights.size(); ++orbital)
          diagnostics::neo_trace("soft_orbital",
              "dimension=%ld,orbital=%ld,block_pairing=%.17g,total=%.17g",
              static_cast<long>(dimension), static_cast<long>(orbital + 1),
              weights[orbital], weights.sum());
      };
    }
    auto problem = std::make_shared<ResponseNeoProblem>(
        reduced_gradient,
        structure_hessian.tangent_size(),
        [&orbital_hessian, &structure_hessian](
            const Eigen::VectorXd& vector) {
          const OrbitalCouplingAction image =
              orbital_hessian.apply_orbital_coupling(vector);
          return ResponseNeoDirection{
              image.orbital_hessian,
              structure_hessian.coordinates(
                  structure_hessian.project(image.scaled_structure_forcing))};
        },
        [&orbital_hessian, &structure_hessian](
            const Eigen::VectorXd& vector) {
          const StructureCouplingAction image =
              structure_hessian.apply_coupling(
                  structure_hessian.expand(vector));
          return ResponseNeoDirection{
              orbital_hessian.apply_structure_coupling_adjoint(
                  image.coefficient_response, image.adjoint_multipliers),
              structure_hessian.coordinates(image.hessian)};
        },
        [&orbital_metric](const Eigen::VectorXd& vector) {
          return orbital_metric.apply(vector);
        },
        [&orbital_metric](const Eigen::VectorXd& covector) {
          return orbital_metric.dual_norm(covector);
        },
        [&chart, &orbital_metric](
            const Eigen::VectorXd& residual,
            double shift,
            double relative_residual_target) {
          if (shift == 0.0) {
            return chart.apply_inverse_reduced_block_preconditioner(residual);
          }
          const auto block_inverse = chart.shifted_block_inverse(shift);
          return apply_shifted_metric_preconditioner(
              residual,
              shift,
              relative_residual_target,
              [&chart](const Eigen::VectorXd& vector) {
                return chart.apply_reduced_curvature(vector);
              },
              [&orbital_metric](const Eigen::VectorXd& vector) {
                return orbital_metric.apply(vector);
              },
              block_inverse);
        },
        std::move(orbital_guess),
        operator_relative_accuracy,
        [&orbital_hessian, &structure_hessian, &chart](
            const Eigen::MatrixXd& vectors) {
          const StructureCoordinateCouplingBlock images =
              structure_hessian.apply_coupling_coordinate_block(vectors);
          ResponseNeoDirectionBlock result;
          result.orbital.resize(
              chart.reduced_size(), vectors.cols());
          result.structure = images.hessian_coordinates;
          for (Eigen::Index column = 0; column < vectors.cols(); ++column) {
            result.orbital.col(column) =
                orbital_hessian.apply_structure_coupling_adjoint(
                    images.coefficient_responses[column],
                    images.adjoint_multipliers[column]);
          }
          return result;
        },
        [&orbital_hessian, &structure_hessian, record,
         operator_relative_accuracy](const Eigen::VectorXd& residual) {
          const StructureResponseBlock response =
              orbital_hessian.solve_structure_response_block(
                  {structure_hessian.expand(residual).scaled_coefficients},
                  operator_relative_accuracy);
          // Count iterative C actions in addition to projected basis actions.
          // Dense response uses the already available complete eigenspectrum.
          record->structure_response_actions += response.block_actions;
          record->coupled_block_actions += response.block_actions;
          return structure_hessian.project_coordinates(
              response.scaled_coefficients.at(0));
        });
    auto workspace = std::make_shared<ResponseNeoWorkspace>(*problem);
    solve_step = [problem = std::move(problem),
                  workspace = std::move(workspace), &observe_candidate,
                  observe_model = std::move(observe_model)](
                     const NeoOptions& neo_options, bool allow_early_trial) {
      // The workspace stores a reference; retain its immutable problem here.
      static_cast<void>(problem);
      return workspace->solve(neo_options, allow_early_trial, observe_candidate, observe_model);
    };
  } else {
    auto problem = std::make_shared<NeoProblem>(
        reduced_gradient,
        [&orbital_hessian](
            const Eigen::VectorXd& vector) {
          return orbital_hessian.apply_orbital(vector);
        },
        [&orbital_metric](const Eigen::VectorXd& vector) {
          return orbital_metric.apply(vector);
        },
        [&orbital_metric](const Eigen::VectorXd& covector) {
          return orbital_metric.dual_norm(covector);
        },
        [&chart](const Eigen::VectorXd& residual) {
          return chart.apply_inverse_reduced_block_preconditioner(residual);
        });
    solve_step = [problem = std::move(problem)](
                     const NeoOptions& neo_options, bool allow_early_trial) {
      const NeoResult orbital = solve_neo(*problem, neo_options, allow_early_trial);
      ResponseNeoResult result;
      result.step = {orbital.step, Eigen::VectorXd{}};
      result.hessian_step = {orbital.hessian_step, Eigen::VectorXd{}};
      result.orbital_metric_step = orbital.metric_step;
      result.kkt_residual = {orbital.kkt_residual, Eigen::VectorXd{}};
      result.minimum_curvature_orbital =
          Eigen::VectorXd::Zero(problem->size());
      result.shift = orbital.shift;
      result.predicted_reduction = orbital.predicted_reduction;
      result.step_norm = orbital.step_norm;
      result.residual_norm = orbital.residual_norm;
      result.residual_target = orbital.residual_target;
      result.curvature_residual_norm = orbital.curvature_residual_norm;
      result.curvature_residual_target = orbital.curvature_residual_target;
      result.iterations = orbital.iterations;
      result.coupled_actions = orbital.hessian_actions;
      result.orbital_actions = orbital.hessian_actions;
      result.boundary = orbital.boundary;
      result.hard_case = orbital.hard_case;
      result.global_curvature_certified =
          orbital.global_curvature_certified;
      result.stop_reason = orbital.stop_reason;
      return result;
    };
  }

  diagnostics::neo_trace_context.keyframe = record->accepted_iteration_index;
  diagnostics::neo_trace_context.attempt = 0;
  // Coordinate-small need not mean physically stationary in an ill-conditioned
  // chart. Route the solve using the same dual norm as the final certificate.
  // TrialReady remains approximate and must pass real energy/gradient checks.
  bool allow_early_trial = neo_trial_allowed(
      physical_gradient_norm, options.gradient_tolerance);
  diagnostics::neo_trace("phase", "physical_gradient=%.17g,coordinate_inf=%.17g,tolerance=%.17g,early_trial=%d",
      physical_gradient_norm, gradient_infinity_norm(reduced_gradient),
      options.gradient_tolerance, static_cast<int>(allow_early_trial));
  while (true) {
    ++diagnostics::neo_trace_context.attempt;
    const double trial_radius = *trust_radius;
    NeoOptions neo_options;
    neo_options.trust_radius = trial_radius;
    // A larger physical norm must not loosen the historical absolute forcing
    // budget for nonlinear trial candidates. Both norms have energy units;
    // the physical norm still defines the certificate. Near stationarity only
    // its invariant budget applies. No acceptance or accuracy threshold changes.
    neo_options.relative_residual_tolerance = neo_forcing_term(
        physical_gradient_norm, chart.reduced_size());
    if (allow_early_trial) {
      const double coordinate_gradient = chart.coordinate_dual_norm(reduced_gradient);
      const double coordinate_budget =
          neo_forcing_term(coordinate_gradient, chart.reduced_size()) *
          coordinate_gradient;
      neo_options.relative_residual_tolerance = std::min(
          neo_options.relative_residual_tolerance,
          coordinate_budget / physical_gradient_norm);
    }
    // Match the outer physical stationarity budget. The solver measures
    // orbital covectors in the dual of this same retraction metric; a small
    // coordinate residual cannot certify an unresolved physical gradient.
    neo_options.absolute_residual_tolerance = options.gradient_tolerance;
    neo_options.maximum_subspace_dimension = 0;
    // Select the lowest-curvature branch of the current microproblem, not
    // merely a stationary Newton solution in a gradient-generated subspace.
    // This is direction construction; outer convergence is unchanged.
    neo_options.require_curvature_certificate = true;
    const ResponseNeoResult step = solve_step(neo_options, allow_early_trial);
    record->micro_iterations += step.iterations;
    record->coupled_block_actions += step.coupled_actions;
    record->orbital_hvp_actions += step.orbital_actions;
    record->structure_response_actions += step.structure_actions;
    if (!step.converged() && !step.trial_ready()) {
      result->termination_reason = "neo_microproblem_not_converged";
      return false;
    }

    const Eigen::VectorXd candidate_parameters =
        build_nonredundant_lifted_trial_parameters(
            accepted_orbitals,
            chart,
            parameter_view,
            step.step.orbital);
    if (is_effectively_zero_step(
            candidate_parameters - current_parameters,
            current_parameters)) {
      result->termination_reason = "neo_zero_orbital_step";
      return false;
    }

    auto trial = objective->evaluate_trial_energy(candidate_parameters, true);
    const double gradient_dot_step =
        reduced_gradient.dot(step.step.orbital);
    const double step_dot_hessian_step =
        step.step.orbital.dot(step.hessian_step.orbital) +
        step.step.structure.dot(step.hessian_step.structure);
    const NeoGlobalizationResult globalization = globalize_neo_trial(
        energy,
        trial.energy,
        gradient_dot_step,
        step_dot_hessian_step,
        step.boundary,
        trial_radius,
        options.minimum_step_size);
    diagnostics::neo_trace("trial",
        "radius=%.17g,next_radius=%.17g,hvp=%d,structure=%d,accepted=%d,rho=%.17g,predicted=%.17g,actual=%.17g",
        trial_radius, globalization.next_radius, step.orbital_actions,
        step.structure_actions, static_cast<int>(globalization.accepted),
        globalization.rho, globalization.predicted_reduction,
        globalization.actual_reduction);
    if (!globalization.accepted) {
      ++record->rejected_trial_count;
      if (!(globalization.next_radius < *trust_radius)) {
        result->termination_reason = "neo_trust_radius_exhausted";
        return false;
      }
      *trust_radius = globalization.next_radius;
      continue;
    }

    if (step.trial_ready()) {
      // Validate a real relaxed gradient before accepting an approximate
      // coupled step. Compare invariant physical norms at their own points:
      // reduced coordinate vectors at different charts are not comparable.
      objective->complete_trial(&trial);
      const double trial_norm = trial_physical_gradient_norm(
          trial, parameter_view, objective->input().complete_active_space);
      const double predicted_norm = orbital_metric.dual_norm(
          reduced_gradient + step.hessian_step.orbital);
      // PySCF's factor-three gradient trust guard, applied to physical norms.
      // This is a magnitude check, not a claimed vector-transport identity.
      const bool reliable = trial_norm <= 3.0 * std::max({
          physical_gradient_norm, predicted_norm, options.gradient_tolerance});
      diagnostics::neo_trace("keyframe", "actual_gradient=%.17g,predicted_gradient=%.17g,reliable=%d",
          trial_norm, predicted_norm, static_cast<int>(reliable));
      if (!reliable) {
        ++record->rejected_trial_count;
        // A rejected gradient keyframe must not discard an energy-model
        // contraction already requested by rho. Never expand a rejected step.
        *trust_radius = std::min(trial_radius, globalization.next_radius);
        diagnostics::neo_trace("keyframe_retry", "radius=%.17g,next_radius=%.17g",
            trial_radius, *trust_radius);
        allow_early_trial = false;
        continue;  // Reuse the retained model at the non-expanded retry radius.
      }
    }

    record->kkt_residual_norm = step.residual_norm;
    record->approximate_step = step.trial_ready();
    record->kkt_residual_target = step.residual_target;
    record->orbital_step_norm = step.step_norm;
    record->shift = step.shift;
    if (options.neo_iteration_callback) {
      record->physical_orbital_kkt_norm =
          orbital_metric.dual_norm(step.kkt_residual.orbital);
    }
    record->curvature_residual_norm = step.curvature_residual_norm;
    record->curvature_residual_target = step.curvature_residual_target;
    record->global_curvature_certified = step.global_curvature_certified;
    record->accepted_trial_radius = trial_radius;
    record->next_trust_radius = globalization.next_radius;
    record->gradient_dot_step += gradient_dot_step;
    record->step_dot_hessian_step += step_dot_hessian_step;
    record->predicted_reduction += globalization.predicted_reduction;
    record->actual_reduction += globalization.actual_reduction;
    record->reached_boundary = record->reached_boundary || step.boundary;
    *trust_radius = globalization.next_radius;

    accepted->parameters =
        parameter_view.pack(trial.orbital_preparation_input);
    accepted->next_orbital_guess_packed =
        chart.expand_step(step.step.orbital);
    accepted->trial = std::move(trial);
    return true;
  }
}

}  // namespace

BackendRunResult run_neo_backend(
    VbScfObjective* objective,
    const SparseParameterLayout& parameter_view,
    const VbScfOptimizerOptions& options,
    const Eigen::VectorXd& initial_parameters,
    const Eigen::VectorXd& initial_gradient,
    double initial_energy,
    VbScfOptimizerResult* result) {
  BackendRunResult run_result;
  Eigen::VectorXd parameters = initial_parameters;
  Eigen::VectorXd gradient = initial_gradient;
  double energy = initial_energy;
  double previous_energy = initial_energy;
  double trust_radius =
      std::max(options.minimum_step_size, options.initial_step_size);
  const OrbitalPreconditioner preconditioner =
      resolve_neo_preconditioner(*objective, options.orbital_preconditioner);
  const bool use_ciah =
      objective->second_order_context()->selected_state_indices.size() == 1;
  OrbitalChart chart = build_orbital_geometry(
      *objective, parameter_view, OrbitalChart::Coordinates::Unwhitened);
  auto projected = chart.project_gradient(gradient);
  double final_gradient_l2 = chart.coordinate_dual_norm(projected.reduced_gradient);
  Eigen::VectorXd orbital_guess_packed;
  Eigen::MatrixXd structure_guess;
  double recycle_gradient_target =
      neo_forcing_term(final_gradient_l2, chart.reduced_size()) *
      final_gradient_l2;

  while (run_result.n_iterations < options.max_iterations) {
    if (chart.reduced_size() == 0) {
      result->termination_reason = "neo_nonredundant_space_empty";
      final_gradient_l2 = 0.0;
      break;
    }
    final_gradient_l2 = chart.coordinate_dual_norm(projected.reduced_gradient);
    if (run_result.n_iterations == 0 &&
        orbital_gradient_converged(*objective, parameter_view, chart,
            projected.reduced_gradient, options.gradient_tolerance)) {
      result->converged = true;
      result->termination_reason = "neo_initial_tolerance";
      break;
    }
    // A terminal point needs geometry and a true gradient, not another model.
    prepare_orbital_curvature(*objective, chart, preconditioner, use_ciah);

    NeoIterationRecord iteration_record;
    iteration_record.accepted_iteration_index = run_result.n_iterations + 1;
    iteration_record.initial_trust_radius = trust_radius;
    AcceptedNeoKeyframe accepted;
    if (!build_accepted_neo_keyframe(
            objective,
            parameter_view,
            options,
            chart,
            projected.reduced_gradient,
            parameters,
            orbital_guess_packed,
            structure_guess,
            energy,
            &trust_radius,
            &iteration_record,
            result,
            &accepted)) {
      break;
    }

    if (!accepted.trial.gradient_ready) objective->complete_trial(&accepted.trial);
    accepted.parameters =
        parameter_view.pack(accepted.trial.orbital_preparation_input);
    gradient = accepted.trial.gradient;
    energy = accepted.trial.energy;
    objective->commit(std::move(accepted.trial));
    parameters = std::move(accepted.parameters);
    orbital_guess_packed = std::move(accepted.next_orbital_guess_packed);
    structure_guess = std::move(accepted.next_structure_guess);
    // CIAH may evaluate several fixed-structure keyframes before this final
    // relaxed point. Include the macro endpoint without discarding that work.
    ++iteration_record.keyframes;

    // Every accepted point has a new exact gradient and a new local model.
    // It is therefore an outer iteration, including for iteration limits and
    // progress output; grouping such points obscures the actual work.
    // Preconditioning changes only curvature surrogates, not the quotient
    // maps. Use one chart for stationarity, reporting, and the next solve.
    OrbitalChart next_geometry = build_orbital_geometry(
        *objective, parameter_view, OrbitalChart::Coordinates::Unwhitened);
    auto next_projected = next_geometry.project_gradient(gradient);
    final_gradient_l2 = next_geometry.coordinate_dual_norm(next_projected.reduced_gradient);
    if (options.neo_iteration_callback) {
      const NonredundantRetractionMetric next_metric(
          next_geometry, parameter_view,
          objective->input().orbital_preparation_input);
      iteration_record.accepted_physical_gradient_norm =
          next_metric.dual_norm(next_projected.reduced_gradient);
    }
    // Retain the packed trial direction only while following the same
    // exact-gradient keyframe sequence. All accepted points are still
    // counted and reported independently.
    if (objective->second_order_context()->selected_state_indices.size() != 1 &&
        (final_gradient_l2 <= recycle_gradient_target ||
         (final_gradient_l2 < options.gradient_tolerance * std::sqrt(next_geometry.reduced_size()) &&
          next_geometry.reference_gradient_inf_norm(next_projected.reduced_gradient) < options.gradient_tolerance))) {
      orbital_guess_packed.resize(0);
      recycle_gradient_target =
          neo_forcing_term(final_gradient_l2, next_geometry.reduced_size()) *
          final_gradient_l2;
    }
    const double energy_change = energy - previous_energy;
    previous_energy = energy;
    if (std::abs(energy_change) < options.energy_tolerance &&
        orbital_gradient_converged(*objective, parameter_view, next_geometry,
            next_projected.reduced_gradient, options.gradient_tolerance)) {
      result->converged = true;
      result->termination_reason = "neo_projected_gradient_tolerance";
    }
    chart = std::move(next_geometry);
    projected = std::move(next_projected);

    ++run_result.n_iterations;
    iteration_record.trust_ratio =
        iteration_record.predicted_reduction > 0.0
            ? iteration_record.actual_reduction /
                  iteration_record.predicted_reduction
            : -std::numeric_limits<double>::infinity();
    result->neo_iteration_trace.push_back(iteration_record);
    if (options.neo_iteration_callback) {
      options.neo_iteration_callback(iteration_record);
    }
    sync_result_from_objective(*objective, result);
    record_accepted_iteration_snapshot(
        objective,
        run_result.n_iterations,
        options,
        result,
        nullptr,
        &projected.reduced_gradient,
        &chart);
    if (result->converged) break;
  }

  result->final_projected_gradient_inf_norm =
      chart.reference_gradient_inf_norm(projected.reduced_gradient);
  result->final_projected_gradient_l2_norm = final_gradient_l2;
  run_result.final_gradient_l2_norm = final_gradient_l2;
  run_result.final_projected_gradient_ready = true;
  return run_result;
}

}  // namespace xmvb::vb::optimizer_detail
