#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#include "input/loading/loader.hpp"
#include "linalg/metric_solve.hpp"
#include "output/trace/binary.hpp"
#include "vbscf/derivatives/gradient/orbital/evaluator.hpp"
#include "vbscf/derivatives/hessian/context/accepted_point.hpp"
#include "vbscf/orbitals/preparation/cholesky.hpp"
#include "vbscf/structures/overlap/representation.hpp"

namespace {

using xmvb::vb::RepresentationOverlap;

void print_spectrum(const char* label, const Eigen::MatrixXd& matrix) {
  const Eigen::MatrixXd symmetric = 0.5 * (matrix + matrix.transpose()).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(
      symmetric, Eigen::EigenvaluesOnly);
  if (eigensolver.info() != Eigen::Success)
    throw std::runtime_error(std::string(label) + " eigensolver failed");
  const double smallest = eigensolver.eigenvalues()[0];
  const double largest = eigensolver.eigenvalues().tail(1)[0];
  std::cout << label << "_lambda_min=" << smallest << '\n'
            << label << "_lambda_max=" << largest << '\n'
            << label << "_condition="
            << (smallest > 0.0 ? largest / smallest
                               : std::numeric_limits<double>::infinity()) << '\n';
}

template<class Action>
void inspect_rhs(const char* label, const Eigen::VectorXd& rhs,
                 const Action& action, const Eigen::VectorXd& jacobi,
                 const Eigen::LLT<Eigen::MatrixXd>& factor) {
  const double rhs_norm = rhs.norm();
  const double target = std::sqrt(std::numeric_limits<double>::epsilon()) * rhs_norm;
  std::cout << label << "_rhs_norm=" << rhs_norm << '\n'
            << label << "_target=" << target << '\n';

  xmvb::linalg::MetricSolveStats stats;
  try {
    const Eigen::VectorXd solution = xmvb::linalg::solve_metric(
        rhs, action,
        [&](const Eigen::VectorXd& residual) {
          return Eigen::VectorXd(residual.cwiseQuotient(jacobi));
        }, &stats);
    const double true_residual = (rhs - action(solution)).norm();
    std::cout << label << "_cg_status=certified\n"
              << label << "_cg_solution_norm=" << solution.norm() << '\n'
              << label << "_cg_true_residual=" << true_residual << '\n'
              << label << "_cg_relative_residual="
              << (rhs_norm > 0.0 ? true_residual / rhs_norm : 0.0) << '\n';
  } catch (const std::exception& error) {
    std::cout << label << "_cg_status=failed\n"
              << label << "_cg_failure=" << error.what() << '\n';
  }
  std::cout << label << "_cg_iterations=" << stats.iterations << '\n'
            << label << "_cg_actions=" << stats.actions << '\n'
            << label << "_cg_certifications=" << stats.certifications << '\n'
            << label << "_cg_last_certified_relative_residual="
            << stats.relative_residual << '\n';

  if (factor.info() != Eigen::Success) {
    std::cout << label << "_llt_status=failed\n";
    return;
  }
  const Eigen::VectorXd solution = factor.solve(rhs);
  if (factor.info() != Eigen::Success || !solution.allFinite()) {
    std::cout << label << "_llt_status=failed\n";
    return;
  }
  const double true_residual = (rhs - action(solution)).norm();
  std::cout << label << "_llt_status="
            << (true_residual <= target ? "certified" : "above_tolerance") << '\n'
            << label << "_llt_solution_norm=" << solution.norm() << '\n'
            << label << "_llt_true_residual=" << true_residual << '\n'
            << label << "_llt_relative_residual="
            << (rhs_norm > 0.0 ? true_residual / rhs_norm : 0.0) << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2 || argc > 3)
      throw std::invalid_argument(
          "usage: test_selected_metric_point INPUT.xmi [OUTPUT_PREFIX]");
    auto loaded = xmvb::vb::load_vbscf_input_with_timings(argv[1]);
    const Eigen::MatrixXd no_initial_eigenvectors;
    xmvb::vb::OrbitalGradientEvaluator evaluator;
    const auto evaluated = evaluator.evaluate_without_reference_energy_gradient(
        loaded.input, {0}, {1.0}, loaded.nuclear_repulsion_energy,
        xmvb::vb::StructureEigensolver::Davidson,
        xmvb::vb::StructureSolveAccuracy{}, no_initial_eigenvectors);
    const auto& context = evaluated.second_order_context;
    if (!context || !context->structure_action)
      throw std::runtime_error("accepted structure action is unavailable");
    const auto& action = *context->structure_action;
    const Eigen::MatrixXd& coefficients = context->selected_state_eigenvectors;
    const auto images = context->selected_state_structure_images
        ? *context->selected_state_structure_images : action.apply(coefficients);
    const auto& orbitals = loaded.input.orbital_preparation_input;
    const int n_active = orbitals.n_active_orbitals;
    const auto& overlap =
        context->prepared_active_space.orbital_result.active_orbital_overlap_matrix;
    const Eigen::MatrixXd upper = xmvb::vb::overlap_upper(
        Eigen::Map<const Eigen::MatrixXd>(overlap.data(), n_active, n_active));
    const auto& alpha =
        context->same_spin_pair_cache.alpha_reuse_table.unique_determinants;
    const auto& beta =
        context->same_spin_pair_cache.beta_reuse_table.unique_determinants;
    const auto prepared = RepresentationOverlap(action, alpha, beta, n_active)
                              .prepare_metric(upper);
    const int dimension = action.n_structures();
    std::cout << std::setprecision(17)
              << "direct_ci=" << action.supports_integral_direction() << '\n'
              << "dimension=" << dimension << '\n'
              << "resident_selected_matrix=" << (prepared.selected_matrix().size() != 0)
              << '\n';
    if (!RepresentationOverlap::PreparedMetric::fits_selected(dimension))
      throw std::runtime_error("selected Gram exceeds bounded diagnostic workspace");

    // The identity block is permitted only by the existing selected-matrix
    // budget. No full spin-carrier or determinant-product matrix is formed.
    const Eigen::MatrixXd metric = prepared.selected_matrix().size()
        ? prepared.selected_matrix()
        : prepared(Eigen::MatrixXd::Identity(dimension, dimension));
    const auto original_action = [&](const Eigen::VectorXd& vector) {
      return Eigen::VectorXd(prepared(vector).col(0));
    };
    const Eigen::MatrixXd symmetric = 0.5 * (metric + metric.transpose()).eval();
    std::cout << "metric_symmetry_relative="
              << (metric - metric.transpose()).norm() / std::max(1.0, metric.norm())
              << '\n';
    print_spectrum("metric", metric);
    const Eigen::VectorXd jacobi = action.preconditioner_diagonal().overlap;
    if (jacobi.size() != dimension || !jacobi.allFinite() ||
        (jacobi.array() <= 0.0).any())
      throw std::runtime_error("accepted Jacobi diagonal is invalid");
    const Eigen::VectorXd inverse_sqrt = jacobi.array().rsqrt();
    const Eigen::MatrixXd scaled =
        inverse_sqrt.asDiagonal() * symmetric * inverse_sqrt.asDiagonal();
    print_spectrum("jacobi_scaled", scaled);
    std::cout << "jacobi_metric_diagonal_relative_gap="
              << (jacobi - metric.diagonal()).norm() /
                     std::max(1.0, metric.diagonal().norm()) << '\n';

    const Eigen::MatrixXd state_h = coefficients.transpose() * images.hamiltonian;
    const Eigen::MatrixXd transport_rhs =
        2.0 * (images.hamiltonian - images.overlap * state_h);
    Eigen::VectorXd deterministic(dimension);
    for (int i = 0; i < dimension; ++i)
      deterministic[i] = std::sin(0.37 * (i + 1)) +
                         std::cos(0.11 * (i + 1));
    const Eigen::VectorXd deterministic_image = original_action(deterministic);
    std::cout << "deterministic_action_matrix_relative_gap="
              << (deterministic_image - metric * deterministic).norm() /
                     std::max(1.0, deterministic_image.norm()) << '\n';

    if (argc == 3) {
      const std::filesystem::path prefix(argv[2]);
      xmvb::output::write_binary_buffer(prefix.string() + ".metric.bin",
                                        metric.data(), metric.size());
      xmvb::output::write_binary_buffer(prefix.string() + ".gc0.bin",
                                        transport_rhs.data(), transport_rhs.size());
      xmvb::output::write_binary_buffer(prefix.string() + ".deterministic.bin",
                                        deterministic.data(), deterministic.size());
    }

    const Eigen::LLT<Eigen::MatrixXd> factor(symmetric);
    std::cout << "llt_factor_success=" << (factor.info() == Eigen::Success) << '\n';
    for (int state = 0; state < transport_rhs.cols(); ++state) {
      const std::string label = "transport_" + std::to_string(state);
      inspect_rhs(label.c_str(), transport_rhs.col(state), original_action,
                  jacobi, factor);
    }
    inspect_rhs("deterministic", deterministic, original_action, jacobi, factor);
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "selected_metric_point_failure=" << error.what() << '\n';
    return 1;
  }
}
