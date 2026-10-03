#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "input/loading/loader.hpp"
#include "vbscf/derivatives/gradient/orbital/evaluator.hpp"
#include "vbscf/derivatives/hessian/exact/operator.hpp"
#include "vbscf/optimization/driver/session.hpp"
#include "vbscf/optimization/objective/function.hpp"
#include "vbscf/optimization/trust_region/retraction.hpp"
#include "vbscf/orbitals/charts/layout.hpp"

namespace {

using namespace xmvb::vb;

double positive_tolerance(const char* text, const char* label) {
  const double value = std::stod(text);
  if (!(value > 0.0) || !std::isfinite(value)) {
    throw std::invalid_argument(std::string(label) + " must be positive and finite");
  }
  return value;
}

Eigen::VectorXd read_packed_direction(
    const std::string& path, Eigen::Index expected_size) {
  std::ifstream stream(path, std::ios::binary | std::ios::ate);
  if (!stream) {
    throw std::runtime_error("cannot open packed direction: " + path);
  }
  const auto expected_bytes =
      static_cast<std::streamoff>(expected_size * sizeof(double));
  if (stream.tellg() != expected_bytes) {
    throw std::runtime_error("packed direction has the wrong byte count: " + path);
  }
  stream.seekg(0);
  Eigen::VectorXd direction(expected_size);
  stream.read(reinterpret_cast<char*>(direction.data()), expected_bytes);
  if (!stream || !direction.allFinite()) {
    throw std::runtime_error("packed direction is incomplete or non-finite: " + path);
  }
  return direction;
}

void write_doubles(const std::string& path, const double* values,
                   std::size_t count) {
  std::ofstream stream(path, std::ios::binary | std::ios::trunc);
  if (!stream) {
    throw std::runtime_error("cannot create diagnostic file: " + path);
  }
  stream.write(reinterpret_cast<const char*>(values),
               static_cast<std::streamsize>(count * sizeof(double)));
  if (!stream) {
    throw std::runtime_error("cannot complete diagnostic file: " + path);
  }
}

double infinity_norm(const Eigen::VectorXd& values) {
  return values.size() == 0 ? 0.0 : values.cwiseAbs().maxCoeff();
}

void run_audit(int argc, char** argv) {
  if (argc < 2 || argc > 7) {
    throw std::invalid_argument(
        "usage: audit_orbital_gradient INPUT [dense|davidson] "
        "[energy-tol] [gradient-tol] [packed-direction-f64.bin] [full-hvp]");
  }
  const bool audit_hvp = argc == 7;
  if (audit_hvp && std::string(argv[6]) != "full-hvp") {
    throw std::invalid_argument("last argument must be full-hvp");
  }

  auto loaded = load_vbscf_input_with_timings(argv[1]);
  StructureEigensolver solver = loaded.structure_eigensolver;
  if (argc >= 3) {
    const std::string name = argv[2];
    if (name == "dense") {
      solver = StructureEigensolver::Dense;
    } else if (name == "davidson") {
      solver = StructureEigensolver::Davidson;
    } else {
      throw std::invalid_argument("eigensolver must be dense or davidson");
    }
  }
  StructureSolveAccuracy accuracy;
  if (argc >= 4) {
    accuracy.energy_tolerance = positive_tolerance(argv[3], "energy tolerance");
  }
  if (argc >= 5) {
    accuracy.gradient_tolerance = positive_tolerance(argv[4], "gradient tolerance");
  }

  if (loaded.state_average_count <= 0) {
    throw std::runtime_error("loaded state average has no selected states");
  }
  std::vector<int> selected_states(loaded.state_average_count);
  std::iota(selected_states.begin(), selected_states.end(), 0);
  const std::vector<double> state_weights(selected_states.size(), 1.0);

  const SparseParameterLayout layout(loaded.input.orbital_preparation_input);
  OrbitalGradientEvaluator evaluator;
  VbScfObjective objective(
      std::move(loaded.input), layout, selected_states, state_weights,
      loaded.nuclear_repulsion_energy, solver, accuracy, &evaluator);
  const Eigen::VectorXd accepted_parameters =
      layout.pack(objective.input().orbital_preparation_input);
  Eigen::VectorXd packed_gradient;
  const double accepted_energy = objective(accepted_parameters, packed_gradient);
  const OrbitalChart chart = optimizer_detail::build_orbital_chart(
      objective, layout, OrbitalPreconditioner::Identity);
  const Eigen::VectorXd reduced_gradient =
      chart.project_reduced_gradient(packed_gradient);
  const NonredundantRetractionMetric metric(
      chart, layout, objective.input().orbital_preparation_input);

  Eigen::VectorXd reduced_direction = -reduced_gradient;
  if (argc >= 6) {
    const Eigen::VectorXd packed_direction =
        read_packed_direction(argv[5], layout.size());
    reduced_direction =
        chart.project_vector(packed_direction).reduced_gradient;
  }
  const double initial_physical_norm = metric.norm(reduced_direction);
  if (!(initial_physical_norm > 0.0) || !std::isfinite(initial_physical_norm)) {
    throw std::runtime_error("probe direction has zero or invalid physical norm");
  }
  reduced_direction /= initial_physical_norm;
  const double analytic_derivative = reduced_gradient.dot(reduced_direction);
  const double raw_derivative =
      packed_gradient.dot(chart.expand_step(reduced_direction));

  const auto& orbital_input = objective.input().orbital_preparation_input;
  const std::string prefix = "audit_orbital_gradient_";
  write_doubles(prefix + "accepted_raw_f64.bin",
                orbital_input.orbital_value_table.data(),
                orbital_input.orbital_value_table.size());
  write_doubles(prefix + "ao_overlap_f64.bin",
                orbital_input.ao_overlap_matrix.data(),
                orbital_input.ao_overlap_matrix.size());
  write_doubles(prefix + "packed_gradient_f64.bin",
                packed_gradient.data(), packed_gradient.size());
  const Eigen::VectorXd packed_tangent = chart.expand_step(reduced_direction);
  write_doubles(prefix + "packed_tangent_f64.bin",
                packed_tangent.data(), packed_tangent.size());

  std::cout << std::setprecision(17)
            << "input = " << argv[1] << '\n'
            << "eigensolver = " << structure_eigensolver_name(solver) << '\n'
            << "two_electron_mode = "
            << standard_two_electron_mode_name(loaded.standard_two_electron_mode) << '\n'
            << "energy_tolerance = " << accuracy.energy_tolerance << '\n'
            << "gradient_tolerance = " << accuracy.gradient_tolerance << '\n'
            << "n_basis_functions = " << orbital_input.n_basis_functions << '\n'
            << "n_orbitals = " << orbital_input.n_orbitals << '\n'
            << "n_active_orbitals = " << orbital_input.n_active_orbitals << '\n'
            << "orbital_type = " << orbital_input.orbital_type << '\n'
            << "complete_active_space = " << objective.input().complete_active_space << '\n'
            << "selected_state_count = " << selected_states.size() << '\n'
            << "accepted_energy = " << accepted_energy << '\n'
            << "packed_parameter_count = " << layout.size() << '\n'
            << "reduced_parameter_count = " << chart.reduced_size() << '\n'
            << "packed_gradient_inf = " << infinity_norm(packed_gradient) << '\n'
            << "reduced_gradient_inf = " << infinity_norm(reduced_gradient) << '\n'
            << "physical_gradient_dual_norm = "
            << metric.dual_norm(reduced_gradient) << '\n'
            << "probe_source = "
            << (argc >= 6 ? argv[5] : "negative_reduced_gradient") << '\n'
            << "probe_physical_norm = " << metric.norm(reduced_direction) << '\n'
            << "analytic_directional_derivative = " << analytic_derivative << '\n'
            << "raw_directional_derivative = " << raw_derivative << '\n'
            << "chart_adjoint_error = "
            << std::abs(raw_derivative - analytic_derivative) << '\n'
            << "binary_layout = float64, Eigen column-major, native byte order\n"
            << "accepted_raw_count = " << orbital_input.orbital_value_table.size() << '\n'
            << "ao_overlap_shape = " << orbital_input.ao_overlap_matrix.rows()
            << " x " << orbital_input.ao_overlap_matrix.cols() << '\n'
            << "packed_gradient_count = " << packed_gradient.size() << '\n'
            << "binary_prefix = " << prefix << '\n';
  std::cout.flush();

  double analytic_curvature = 0.0;
  if (audit_hvp) {
    ExactHvpOperator operation(objective.second_order_context(),
                              &objective.input(), layout, &chart);
    const Eigen::VectorXd image = operation.apply_reduced(reduced_direction);
    analytic_curvature = reduced_direction.dot(image);
    HvpComponents fixed_ci;
    fixed_ci.structure_response = false;
    const Eigen::VectorXd fixed_image =
        operation.apply_reduced(reduced_direction, fixed_ci);
    const auto diagnostics = operation.diagnostics();
    std::cout << "full_hvp_curvature = " << analytic_curvature << '\n'
              << "fixed_ci_hvp_curvature = "
              << reduced_direction.dot(fixed_image) << '\n'
              << "structure_response_curvature = "
              << reduced_direction.dot(image - fixed_image) << '\n'
              << "hvp_response_enabled = " << diagnostics.outer_response_enabled << '\n'
              << "hvp_response_relative_residual = "
              << diagnostics.max_structure_response_relative_residual << '\n';
    write_doubles(prefix + "reduced_hvp_f64.bin", image.data(), image.size());
    std::cout.flush();
  }
  const std::vector<double> steps = audit_hvp
      ? std::vector<double>{1.0e-2, 3.0e-3, 1.0e-3, 3.0e-4}
      : std::vector<double>{1.0e-3, 3.0e-4};
  for (const double h : steps) {
    const auto energy_at = [&](double signed_step) {
      const OrbitalPreparationInput trial_orbitals = chart.retract_step(
          orbital_input, reduced_direction, signed_step);
      const Eigen::VectorXd trial_parameters = layout.pack(trial_orbitals);
      return objective.evaluate_trial_energy(trial_parameters, false).energy;
    };
    const double plus_energy = energy_at(h);
    const double minus_energy = energy_at(-h);
    const double finite_difference = (plus_energy - minus_energy) / (2.0 * h);
    const double error = finite_difference - analytic_derivative;
    std::cout << "probe_h = " << h << '\n'
              << "probe_plus_energy = " << plus_energy << '\n'
              << "probe_minus_energy = " << minus_energy << '\n'
              << "probe_fd_derivative = " << finite_difference << '\n'
              << "probe_fd_minus_analytic = " << error << '\n'
              << "probe_fd_relative_error = "
              << std::abs(error) /
                     std::max({1.0e-12, std::abs(finite_difference),
                               std::abs(analytic_derivative)})
              << '\n';
    if (audit_hvp) {
      const double curvature =
          ((plus_energy - accepted_energy) +
           (minus_energy - accepted_energy)) / (h * h);
      std::cout << "probe_fd_curvature = " << curvature << '\n'
                << "probe_curvature_minus_hvp = "
                << curvature - analytic_curvature << '\n';
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    run_audit(argc, argv);
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "audit_orbital_gradient: " << error.what() << '\n';
    return 1;
  }
}
