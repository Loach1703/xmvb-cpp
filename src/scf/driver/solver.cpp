#include "scf/driver/solver.hpp"

#include <cmath>
#include <span>
#include <stdexcept>
#include <utility>

#include "linalg/eigensolver.hpp"
#include "scf/mixing/diis.hpp"

namespace xmvb::scf {

namespace {

Eigen::MatrixXd symmetrize(
    const Eigen::Ref<const Eigen::MatrixXd>& matrix) {
  return 0.5 * (matrix + matrix.transpose()).eval();
}

void diagonalize(
    const linalg::GeneralizedEigensolver& eigensolver,
    const Eigen::Ref<const Eigen::MatrixXd>& fock,
    const Eigen::Ref<const Eigen::MatrixXd>& overlap,
    Result* result) {
  const int n_basis = static_cast<int>(fock.rows());
  auto eigenpairs = eigensolver.solve_dense(
      std::span<const double>(fock.data(), fock.size()),
      std::span<const double>(overlap.data(), overlap.size()),
      n_basis);
  result->orbital_energies = Eigen::Map<const Eigen::VectorXd>(
      eigenpairs.eigenvalues.data(), n_basis);
  result->orbitals = Eigen::Map<const Eigen::MatrixXd>(
      eigenpairs.eigenvector_matrix.data(), n_basis, n_basis);
}

Eigen::MatrixXd occupied_density(
    const Eigen::Ref<const Eigen::MatrixXd>& orbitals,
    int n_occupied) {
  return 2.0 * (orbitals.leftCols(n_occupied) *
                orbitals.leftCols(n_occupied).transpose()).eval();
}

double max_abs_difference(
    const Eigen::Ref<const Eigen::MatrixXd>& left,
    const Eigen::Ref<const Eigen::MatrixXd>& right) {
  return (left - right).cwiseAbs().maxCoeff();
}

Evaluation checked_evaluation(
    const Evaluate& evaluate,
    const Eigen::Ref<const Eigen::MatrixXd>& density,
    int n_basis) {
  Evaluation physical = evaluate(density);
  if (physical.fock.rows() != n_basis || physical.fock.cols() != n_basis) {
    throw std::invalid_argument("SCF evaluation returned an invalid Fock shape");
  }
  if (!physical.fock.allFinite() ||
      !std::isfinite(physical.electronic_energy)) {
    throw std::runtime_error("SCF evaluation returned a non-finite result");
  }
  return physical;
}

}  // namespace

Result solve_restricted(
    int n_electrons,
    const Eigen::Ref<const Eigen::MatrixXd>& overlap,
    const Eigen::Ref<const Eigen::MatrixXd>& core,
    const Evaluate& evaluate,
    const Options& options,
    const Eigen::MatrixXd& initial_density) {
  const int n_basis = static_cast<int>(core.rows());
  if (n_basis <= 0 || core.cols() != n_basis ||
      overlap.rows() != n_basis || overlap.cols() != n_basis) {
    throw std::invalid_argument("SCF AO matrix dimensions do not match");
  }
  if (n_electrons < 0 || n_electrons % 2 != 0 ||
      n_electrons > 2 * n_basis) {
    throw std::invalid_argument("restricted SCF requires an even electron count within basis capacity");
  }
  if (!evaluate) {
    throw std::invalid_argument("SCF evaluation callback must be set");
  }
  if (options.max_iterations <= 0 || !(options.density_tolerance > 0.0) ||
      options.old_density_weight < 0.0 ||
      options.old_density_weight >= 1.0 ||
      options.diis_history < 0 || options.diis_start_iteration < 0) {
    throw std::invalid_argument("invalid restricted SCF options");
  }
  if (initial_density.size() != 0 &&
      (initial_density.rows() != n_basis || initial_density.cols() != n_basis)) {
    throw std::invalid_argument("initial SCF density dimensions do not match");
  }

  const Eigen::MatrixXd symmetric_overlap = symmetrize(overlap);
  const Eigen::MatrixXd symmetric_core = symmetrize(core);
  linalg::GeneralizedEigensolver eigensolver;
  Result result;
  diagonalize(eigensolver, symmetric_core, symmetric_overlap, &result);
  Eigen::MatrixXd density = initial_density.size() != 0
      ? symmetrize(initial_density).eval()
      : occupied_density(result.orbitals, n_electrons / 2);
  Diis diis(options.diis_history);

  for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
    Evaluation physical = checked_evaluation(evaluate, density, n_basis);
    const Eigen::MatrixXd error =
        commutator_error(physical.fock, density, symmetric_overlap);
    diis.push(physical.fock, error);

    Eigen::MatrixXd orbital_fock;
    const bool use_diis =
        iteration + 1 >= options.diis_start_iteration &&
        diis.extrapolate(&orbital_fock);
    if (!use_diis) {
      orbital_fock = physical.fock;
    }
    diagonalize(eigensolver, orbital_fock, symmetric_overlap, &result);
    Eigen::MatrixXd candidate =
        occupied_density(result.orbitals, n_electrons / 2);
    result.iterations = iteration + 1;

    const bool density_close =
        max_abs_difference(candidate, density) < options.density_tolerance;
    const bool residual_close =
        error.cwiseAbs().maxCoeff() < options.density_tolerance;
    const bool last_iteration = iteration + 1 == options.max_iterations;
    if ((density_close && residual_close) || last_iteration) {
      // Report a physical evaluation at the returned orbital density. DIIS
      // changes only the diagonalized Fock, never the physical energy.
      Evaluation candidate_physical =
          checked_evaluation(evaluate, candidate, n_basis);
      const double candidate_residual = commutator_error(
          candidate_physical.fock, candidate, symmetric_overlap)
              .cwiseAbs().maxCoeff();
      result.converged = density_close && residual_close &&
          candidate_residual < options.density_tolerance;
      if (result.converged || last_iteration) {
        result.density = std::move(candidate);
        result.fock = std::move(candidate_physical.fock);
        result.electronic_energy = candidate_physical.electronic_energy;
        return result;
      }
    }

    density = options.old_density_weight * density +
        (1.0 - options.old_density_weight) * candidate;
  }

  throw std::logic_error("restricted SCF iteration loop exited unexpectedly");
}

}  // namespace xmvb::scf
