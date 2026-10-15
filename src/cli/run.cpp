#include "cli/run.hpp"

#include <chrono>
#include <filesystem>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "output/text/report.hpp"
#include "output/trace/binary.hpp"
#include "output/molden/writer.hpp"
#include "output/trace/accepted_iteration.hpp"
#include "output/trace/neo.hpp"
#include "output/trace/tnhvp.hpp"
#include "input/loading/loader.hpp"
#include "vbscf/optimization/driver/optimizer.hpp"
#include "input/deck/keywords.hpp"
#ifdef XMVB_ENABLE_TNVB
#include "tnvb/electronic/evaluator.hpp"
#endif

namespace xmvb::cli {

namespace fs = std::filesystem;

int run(Options command_line) {
  const std::string& input_path = command_line.input_path;
  auto& load_options = command_line.load;
  auto& options = command_line.optimizer;
  const std::string& dump_trace_dir = command_line.trace_directory;
  const std::string& tnhvp_trace_path = command_line.tnhvp_trace_path;
  const std::string& neo_trace_path = command_line.neo_trace_path;
  const std::string& dump_final_orbital_value_table_bin =
      command_line.final_orbitals_path;
  const bool user_specified_max_iterations =
      command_line.max_iterations_explicit;

  const auto command_start_time = std::chrono::system_clock::now();
  const auto command_start_steady_time = std::chrono::steady_clock::now();
  auto load_result =
      xmvb::vb::load_vbscf_input_with_timings(input_path, load_options);
  if (!command_line.optimizer_backend_explicit) {
    switch (load_result.scf_optimizer) {
      case xmvb::vb::InputScfOptimizer::Lbfgs:
        options.backend = xmvb::vb::VbScfOptimizerBackend::Lbfgs;
        break;
      case xmvb::vb::InputScfOptimizer::Tnhvp:
        options.backend =
            xmvb::vb::VbScfOptimizerBackend::NonredundantTruncatedNewton;
        break;
      case xmvb::vb::InputScfOptimizer::BlockLbfgs:
        options.backend = xmvb::vb::VbScfOptimizerBackend::BlockLbfgs;
        break;
      case xmvb::vb::InputScfOptimizer::Neo:
        options.backend = xmvb::vb::VbScfOptimizerBackend::Neo;
        break;
      case xmvb::vb::InputScfOptimizer::Unspecified:
        break;
    }
  }
  if (!command_line.structure_eigensolver_explicit) {
    options.structure_eigensolver = load_result.structure_eigensolver;
  }
  if (!user_specified_max_iterations) {
    // Keep the standalone SCF loop aligned with the input deck semantics:
    // `.xmi` `itmax` controls the maximum iteration count, and omitted `itmax`
    // falls back to the project default of 2000.
    options.max_iterations = load_result.requested_scf_max_iterations;
  }
  const bool use_tnvb = load_result.input.tnvb.enabled;
  std::shared_ptr<const xmvb::vb::ElectronicStateBackend> electronic_backend;
#ifdef XMVB_ENABLE_TNVB
  std::shared_ptr<xmvb::tnvb::ElectronicEvaluator> tnvb_backend;
#endif
  if (use_tnvb) {
#ifndef XMVB_ENABLE_TNVB
    throw std::invalid_argument("this build does not provide the TNVB backend (XMVB_ENABLE_TNVB=OFF)");
#else
    const auto& orbitals=load_result.input.orbital_preparation_input;
    if (options.max_iterations < 0 || load_result.state_average_count != 1 ||
        (options.backend != xmvb::vb::VbScfOptimizerBackend::Lbfgs &&
         options.backend != xmvb::vb::VbScfOptimizerBackend::BlockLbfgs) ||
        options.structure_eigensolver != xmvb::vb::StructureEigensolver::Dense ||
        (orbitals.orbital_type != xmvb::vb::kOrbitalTypeHao &&
         orbitals.orbital_type != xmvb::vb::kOrbitalTypeOeo)) {
      throw std::invalid_argument("TNVB supports ISCF=5/8, one singlet state, HAO/OEO and dense H/S");
    }
    if (options.max_iterations > 0 && load_result.input.tnvb.total_bond <
        load_result.input.tnvb.structure_bond) {
      throw std::invalid_argument("TNVB optimization requires explicit TNSUMBOND >= TNBOND");
    }
    tnvb_backend=std::make_shared<xmvb::tnvb::ElectronicEvaluator>();
    electronic_backend=tnvb_backend;
#endif
  }

  xmvb::output::print_header(
      input_path,
      options,
      load_result,
      command_start_time);

  std::shared_ptr<xmvb::output::AcceptedIterationTraceWriter> trace_writer;
  if (options.verbose && options.max_iterations > 0) {
    options.accepted_iteration_callback = xmvb::output::combine_callbacks(
        std::move(options.accepted_iteration_callback),
        xmvb::output::iteration_logger(options.backend));
  }
  if (!dump_trace_dir.empty()) {
    trace_writer = std::make_shared<xmvb::output::AcceptedIterationTraceWriter>(
        dump_trace_dir,
        input_path,
        load_result,
        xmvb::vb::vbscf_optimizer_backend_name(options.backend));
    options.retain_accepted_iteration_trace = false;
    options.accepted_iteration_callback_requires_reference_gradient = true;
    options.accepted_iteration_callback_requires_full_snapshot = true;
    options.accepted_iteration_callback = xmvb::output::combine_callbacks(
        std::move(options.accepted_iteration_callback),
        [trace_writer](const xmvb::vb::VbScfAcceptedIterationSnapshot& snapshot) {
          trace_writer->write_accepted_iteration(snapshot);
        });
  }

  xmvb::vb::VbScfOptimizer optimizer(options,electronic_backend);
  std::vector<int> selected_states(load_result.state_average_count);
  std::iota(selected_states.begin(), selected_states.end(), 0);
  const std::vector<double> equal_weights(
      selected_states.size(),
      1.0);
  xmvb::vb::VbScfOptimizerResult result;
  if (use_tnvb && options.max_iterations == 0) {
#ifdef XMVB_ENABLE_TNVB
    xmvb::vb::ElectronicStateRequest request;
    request.nuclear_repulsion_energy=load_result.nuclear_repulsion_energy;
    request.eigensolver=xmvb::vb::StructureEigensolver::Dense;
    auto evaluation=tnvb_backend->evaluate_energy_result(load_result.input,request);
    result.energy_only=true;
    result.termination_reason="energy_only";
    result.scf_result=std::move(evaluation.scf_result);
    result.initial_total_energy=result.final_total_energy=result.scf_result.total_energy;
    result.initial_one_electron_reference_energy=result.final_one_electron_reference_energy=
        result.scf_result.one_electron_reference_energy;
    result.total_energy_history={result.final_total_energy};
    result.total_wall_time_seconds=evaluation.total_wall_time_seconds;
    result.optimized_input=std::move(load_result.input);
#endif
  } else {
    result=optimizer.optimize(
          std::move(load_result.input),
          selected_states,
          equal_weights,
          load_result.nuclear_repulsion_energy);
  }
  if (!tnhvp_trace_path.empty()) {
    xmvb::output::write_tnhvp_trace(tnhvp_trace_path, result);
  }
  if (!neo_trace_path.empty()) {
    xmvb::output::write_neo_trace(neo_trace_path, result);
  }
  if (trace_writer != nullptr) {
    trace_writer->finalize(result);
  }
  if (!dump_final_orbital_value_table_bin.empty()) {
    // The final orbital table is the minimal state needed for reduced-chart
    // finite-difference diagnostics.  Keep this separate from trace dumping so
    // production convergence tests do not pay for per-iteration matrix dumps.
    xmvb::output::write_binary_container(
        fs::path(dump_final_orbital_value_table_bin),
        result.optimized_input.orbital_preparation_input.orbital_value_table);
  }
  const bool command_converged = result.converged || result.energy_only;
  std::optional<fs::path> molden_output_path;
  if (load_result.request_molden_output) {
    molden_output_path =
        xmvb::vb::write_molden_file(fs::path(input_path), result.optimized_input);
  }

  std::optional<fs::path> trace_sample_directory;
  if (trace_writer != nullptr) {
    trace_sample_directory = trace_writer->sample_directory();
  }
  xmvb::output::print_summary(
      options,
      load_result,
      result,
      trace_sample_directory,
      molden_output_path,
      command_start_time,
      command_start_steady_time);

  return command_converged ? 0 : 2;
}

}  // namespace xmvb::cli
