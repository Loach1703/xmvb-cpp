#include <algorithm>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "input/loading/loader.hpp"
#include "vbscf/core/storage/eigen.hpp"
#include "vbscf/derivatives/gradient/active_space/evaluator.hpp"
#include "vbscf/integrals/active/two_electron/response/backpropagator.hpp"
#include "vbscf/integrals/active/two_electron/response/adjoint.hpp"
#include "vbscf/integrals/active/two_electron/response/directional.hpp"
#include "vbscf/integrals/ao/one_electron/backpropagator.hpp"

namespace {

struct Options {
  std::string input_path;
  double step = 1.0e-6;
};

void print_usage() {
  std::cerr
      << "usage: check_exact_two_electron_hvp <input.xmi> [--step h]\n";
}

Options parse_arguments(int argc, char** argv) {
  if (argc < 2 || ((argc - 2) % 2 != 0)) {
    print_usage();
    throw std::invalid_argument("invalid argument count");
  }

  Options options;
  options.input_path = argv[1];
  for (int argument_index = 2; argument_index < argc; argument_index += 2) {
    const std::string name = argv[argument_index];
    const std::string value = argv[argument_index + 1];
    if (name == "--step") {
      options.step = std::stod(value);
      continue;
    }
    throw std::invalid_argument("unknown argument: " + name);
  }
  if (!(options.step > 0.0)) {
    throw std::invalid_argument("--step must be positive");
  }
  return options;
}

Eigen::MatrixXd build_dense_active_direction(
    int n_basis_functions,
    int n_active_orbitals) {
  Eigen::MatrixXd direction =
      Eigen::MatrixXd::Zero(n_basis_functions, n_active_orbitals);
  double norm_sq = 0.0;
  for (int active_orbital_index = 0;
       active_orbital_index < n_active_orbitals;
       ++active_orbital_index) {
    for (int basis_function_index = 0;
         basis_function_index < n_basis_functions;
         ++basis_function_index) {
      const std::size_t linear_index =
          active_orbital_index * n_basis_functions +
          basis_function_index;
      const double value =
          std::sin(0.37 * static_cast<double>(linear_index + 1)) +
          0.5 * std::cos(0.13 * static_cast<double>(linear_index + 3));
      direction(basis_function_index, active_orbital_index) = value;
      norm_sq += value * value;
    }
  }
  const double norm = std::sqrt(norm_sq);
  direction /= norm;
  return direction;
}

double max_abs_difference(
    const Eigen::Ref<const Eigen::MatrixXd>& left,
    const Eigen::Ref<const Eigen::MatrixXd>& right) {
  if (left.rows() != right.rows() || left.cols() != right.cols()) {
    throw std::invalid_argument("matrix shape mismatch");
  }
  double max_abs_diff = 0.0;
  for (int col = 0; col < left.cols(); ++col) {
    for (int row = 0; row < left.rows(); ++row) {
      max_abs_diff = std::max(
          max_abs_diff,
          std::abs(left(row, col) - right(row, col)));
    }
  }
  return max_abs_diff;
}

double max_abs_value(const Eigen::Ref<const Eigen::MatrixXd>& values) {
  double max_abs = 0.0;
  for (int col = 0; col < values.cols(); ++col) {
    for (int row = 0; row < values.rows(); ++row) {
      max_abs = std::max(max_abs, std::abs(values(row, col)));
    }
  }
  return max_abs;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse_arguments(argc, argv);

    xmvb::vb::VbScfInputLoadOptions load_options;
    load_options.standard_two_electron_mode =
        xmvb::vb::StandardTwoElectronMode::Exact;
    const auto load_result =
        xmvb::vb::load_vbscf_input_with_timings(options.input_path, load_options);
    const auto& input = load_result.input;

    xmvb::vb::ActiveSpaceGradientEvaluator evaluator;
    const auto gradient_result =
        evaluator.evaluate(
            input,
            {0},
            {1.0},
            load_result.nuclear_repulsion_energy);

    const int n_basis_functions = input.orbital_preparation_input.n_basis_functions;
    const int n_active_orbitals = input.orbital_preparation_input.n_active_orbitals;
    const int n_inactive_doubly_occupied_orbitals =
        (input.orbital_preparation_input.n_total_electrons -
         input.orbital_preparation_input.n_active_electrons) /
        2;
    if (n_basis_functions <= 0 || n_active_orbitals <= 0) {
      throw std::runtime_error("active-space dimensions must be positive");
    }

    const auto dense_active_direction =
        build_dense_active_direction(
            n_basis_functions,
            n_active_orbitals);
    auto accepted_cache = xmvb::vb::build_exact_packed_active_two_electron_adjoint_cache(
        gradient_result.packed_active_two_electron_gradient,
        gradient_result.active_space_two_electron_result.dense_active_coefficients,
        input.ao_integral_input, n_active_orbitals,
        gradient_result.active_space_two_electron_result);
    const auto contracted_gradient = xmvb::vb::backpropagate_exact_packed_active_two_electron_gradient(
        gradient_result.packed_active_two_electron_gradient, accepted_cache);
    const auto contracted_direction =
        xmvb::vb::compute_exact_packed_active_two_electron_integral_directional_derivative_batch(
            accepted_cache, {dense_active_direction}, input.ao_integral_input);
    const auto measure = [&](const char* name) {
      constexpr int repeats = 20;
      const auto start = std::chrono::steady_clock::now();
      double checksum = 0.0;
      for (int repeat = 0; repeat < repeats; ++repeat) {
        const auto gradient = xmvb::vb::backpropagate_exact_packed_active_two_electron_gradient(
            gradient_result.packed_active_two_electron_gradient, accepted_cache);
        const auto integrals =
            xmvb::vb::compute_exact_packed_active_two_electron_integral_directional_derivative_batch(
                accepted_cache, {dense_active_direction}, input.ao_integral_input);
        checksum += gradient.squaredNorm() + integrals.squaredNorm();
      }
      std::cout << name << "_seconds_per_call = "
                << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() / repeats
                << "\n" << name << "_checksum = " << checksum << '\n';
    };
    measure("contracted_kernel");
    // Independent AO-pair contraction reference; move, never duplicate, the image.
    auto orbital_products = std::move(accepted_cache.orbital_pair_products);
    measure("ao_pair_kernel");
    const auto pair_gradient = xmvb::vb::backpropagate_exact_packed_active_two_electron_gradient(
        gradient_result.packed_active_two_electron_gradient, accepted_cache);
    const auto pair_direction =
        xmvb::vb::compute_exact_packed_active_two_electron_integral_directional_derivative_batch(
            accepted_cache, {dense_active_direction}, input.ao_integral_input);
    const double contracted_gradient_error = max_abs_difference(contracted_gradient, pair_gradient);
    const double contracted_direction_error = max_abs_difference(contracted_direction, pair_direction);
    if (contracted_gradient_error > 1.e-11 * std::max(1., max_abs_value(pair_gradient)) ||
        contracted_direction_error > 1.e-11 * std::max(1., max_abs_value(pair_direction))) {
      throw std::runtime_error("one-orbital integral contractions disagree with AO-pair reference");
    }
    const auto analytic_dense_active_gradient_direction =
        xmvb::vb::apply_exact_packed_active_two_electron_adjoint_hessian_vector(
            gradient_result.packed_active_two_electron_gradient,
            gradient_result.active_space_two_electron_result
                .dense_active_coefficients,
            dense_active_direction,
            input.ao_integral_input,
            n_active_orbitals,
            gradient_result.active_space_two_electron_result);
    auto exchange_image = std::move(accepted_cache.exchange_integrals);
    const auto pair_hvp = xmvb::vb::apply_exact_packed_active_two_electron_adjoint_hessian_vector(
        accepted_cache, dense_active_direction, input.ao_integral_input);
    accepted_cache.exchange_integrals = std::move(exchange_image);
    const auto exchange_hvp = xmvb::vb::apply_exact_packed_active_two_electron_adjoint_hessian_vector(
        accepted_cache, dense_active_direction, input.ao_integral_input);
    const double exchange_error = max_abs_difference(exchange_hvp, pair_hvp);
    std::cout << "exchange_integral_bytes = " << accepted_cache.exchange_integrals.size() * sizeof(double)
              << "\nexchange_hvp_reference_error = " << exchange_error << '\n';
    if (exchange_error > 1.0e-11 * std::max(1.0, max_abs_value(pair_hvp)))
      throw std::runtime_error("exchange HVP disagrees with AO-pair reference");
    if (accepted_cache.exchange_integrals.size() != 0) {
      Eigen::MatrixXd active(n_active_orbitals, n_active_orbitals);
      for (int a = 0; a < n_active_orbitals; ++a)
        for (int b = 0; b < n_active_orbitals; ++b)
          active(a, b) = std::sin(0.17 * (a + 1) + 0.31 * (b + 1));
      const auto& coefficients = *accepted_cache.accepted_active_coefficients;
      const Eigen::MatrixXd ao_gradient = coefficients * active.transpose() * coefficients.transpose();
      const auto reference = xmvb::vb::AoEffectiveOneElectronBackpropagator().backpropagate(
          ao_gradient, input.ao_integral_input).inactive_density_gradient;
      const auto contracted = xmvb::vb::backpropagate_exact_active_h1e_gradient(
          xmvb::vb::flatten_matrix_column_major(active), accepted_cache);
      const int n_bf = accepted_cache.n_basis_functions;
      const Eigen::Map<const Eigen::MatrixXd> r(reference.data(), n_bf, n_bf);
      const Eigen::Map<const Eigen::MatrixXd> c(contracted.data(), n_bf, n_bf);
      const Eigen::MatrixXd physical_reference = r + r.transpose();
      const Eigen::MatrixXd physical_contracted = c + c.transpose();
      const double error = max_abs_difference(physical_reference, physical_contracted);
      std::cout << "active_h1e_adjoint_reference_error = " << error << '\n';
      if (error > 1.0e-11 * std::max(1.0, max_abs_value(physical_reference)))
        throw std::runtime_error("active H1e adjoint disagrees with AO-ERI transpose");
    }
    const auto measure_hvp = [&](const char* label) {
      xmvb::vb::ExactPackedActiveTwoElectronApplyWorkspace workspace;
      Eigen::MatrixXd result;
      constexpr int repeats = 10;
      double checksum = 0.0;
      const auto start = std::chrono::steady_clock::now();
      for (int repeat = 0; repeat < repeats; ++repeat) {
        xmvb::vb::apply_exact_packed_active_two_electron_adjoint_hessian_vector(
            accepted_cache, dense_active_direction, input.ao_integral_input,
            &workspace, &result);
        checksum += result.squaredNorm();
      }
      std::cout << label << "_seconds_per_call = "
                << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() / repeats
                << "\n" << label << "_checksum = " << checksum << '\n';
    };
    measure_hvp("exchange_hvp");
    exchange_image = std::move(accepted_cache.exchange_integrals);
    measure_hvp("ao_pair_hvp");
    accepted_cache.exchange_integrals = std::move(exchange_image);
    auto physical_adjoint = std::move(accepted_cache.active_pair_gradient_matrix);
    const Eigen::Index pairs = physical_adjoint.rows();
    accepted_cache.active_pair_gradient_matrix.resize(pairs, pairs);
    for (Eigen::Index p = 0; p < pairs; ++p) {
      for (Eigen::Index q = 0; q <= p; ++q) {
        const double value = std::sin(0.17 * (p + 1) + 0.23 * (q + 1));
        accepted_cache.active_pair_gradient_matrix(p, q) = value;
        accepted_cache.active_pair_gradient_matrix(q, p) = value;
      }
    }
    const auto arbitrary_exchange = xmvb::vb::apply_exact_packed_active_two_electron_adjoint_hessian_vector(
        accepted_cache, dense_active_direction, input.ao_integral_input);
    exchange_image = std::move(accepted_cache.exchange_integrals);
    const auto arbitrary_pair = xmvb::vb::apply_exact_packed_active_two_electron_adjoint_hessian_vector(
        accepted_cache, dense_active_direction, input.ao_integral_input);
    accepted_cache.exchange_integrals = std::move(exchange_image);
    accepted_cache.active_pair_gradient_matrix = std::move(physical_adjoint);
    const double arbitrary_error = max_abs_difference(arbitrary_exchange, arbitrary_pair);
    std::cout << "arbitrary_exchange_hvp_error = " << arbitrary_error << '\n';
    if (arbitrary_error > 1.0e-11 * std::max(1.0, max_abs_value(arbitrary_pair)))
      throw std::runtime_error("exchange HVP failed arbitrary-adjoint reference");
    auto streamed_active_space_result =
        gradient_result.active_space_two_electron_result;
    streamed_active_space_result.dense_ao_pair_products.resize(0, 0);
    const auto streamed_dense_active_gradient_direction =
        xmvb::vb::apply_exact_packed_active_two_electron_adjoint_hessian_vector(
            gradient_result.packed_active_two_electron_gradient,
            streamed_active_space_result.dense_active_coefficients,
            dense_active_direction,
            input.ao_integral_input,
            n_active_orbitals,
            streamed_active_space_result);
    const double streamed_max_abs_diff =
        max_abs_difference(
            analytic_dense_active_gradient_direction,
            streamed_dense_active_gradient_direction);
    if (streamed_max_abs_diff > 1.0e-11) {
      throw std::runtime_error(
          "streamed and resident exact 2e HVPs do not agree");
    }

    std::vector<double> plus_auxiliary_matrix =
        xmvb::vb::flatten_matrix_column_major(
            gradient_result.orbital_preparation_result.auxiliary_orbital_matrix);
    std::vector<double> minus_auxiliary_matrix =
        plus_auxiliary_matrix;
    for (int basis_function_index = 0;
         basis_function_index < n_basis_functions;
         ++basis_function_index) {
      for (int active_orbital_index = 0;
           active_orbital_index < n_active_orbitals;
           ++active_orbital_index) {
        const int column_index =
            n_inactive_doubly_occupied_orbitals + active_orbital_index;
        const std::size_t flat_index =
            column_index * n_basis_functions + basis_function_index;
        const double delta =
            options.step *
            dense_active_direction(basis_function_index, active_orbital_index);
        plus_auxiliary_matrix[flat_index] += delta;
        minus_auxiliary_matrix[flat_index] -= delta;
      }
    }

    xmvb::vb::ActiveSpaceTwoElectronBackpropagator backpropagator;
    const auto plus_backpropagation_result =
        backpropagator.backpropagate(
            gradient_result.packed_active_two_electron_gradient,
            input.ao_integral_input.pair_graph,
            plus_auxiliary_matrix,
            n_basis_functions,
            n_inactive_doubly_occupied_orbitals,
            n_active_orbitals);
    const auto minus_backpropagation_result =
        backpropagator.backpropagate(
            gradient_result.packed_active_two_electron_gradient,
            input.ao_integral_input.pair_graph,
            minus_auxiliary_matrix,
            n_basis_functions,
            n_inactive_doubly_occupied_orbitals,
            n_active_orbitals);

    const auto& plus_dense_active_gradient =
        plus_backpropagation_result.active_auxiliary_orbital_gradient;
    const auto& minus_dense_active_gradient =
        minus_backpropagation_result.active_auxiliary_orbital_gradient;
    const Eigen::MatrixXd finite_difference_dense_active_gradient_direction =
        (plus_dense_active_gradient - minus_dense_active_gradient) /
        (2.0 * options.step);

    const double max_abs_diff =
        max_abs_difference(
            analytic_dense_active_gradient_direction,
            finite_difference_dense_active_gradient_direction);
    const double max_abs_fd =
        max_abs_value(finite_difference_dense_active_gradient_direction);
    std::cout << std::setprecision(12);
    std::cout << "input = " << options.input_path << '\n';
    std::cout << "n_basis_functions = " << n_basis_functions << '\n';
    std::cout << "n_active_orbitals = " << n_active_orbitals << '\n';
    std::cout << "finite_difference_step = " << options.step << '\n';
    std::cout << "analytic_inf_norm = "
              << max_abs_value(analytic_dense_active_gradient_direction) << '\n';
    std::cout << "fd_inf_norm = " << max_abs_fd << '\n';
    std::cout << "max_abs_diff = " << max_abs_diff << '\n';
    std::cout << "max_rel_diff = "
              << max_abs_diff / std::max(1.0, max_abs_fd) << '\n';
    std::cout << "streamed_max_abs_diff = " << streamed_max_abs_diff << '\n';
    std::cout << "contracted_gradient_max_abs_diff = " << contracted_gradient_error << '\n';
    std::cout << "contracted_direction_max_abs_diff = " << contracted_direction_error << '\n';
    std::cout << "orbital_pair_products_bytes = " << orbital_products.size() * sizeof(double) << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
