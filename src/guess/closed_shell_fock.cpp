#include "guess/closed_shell_fock.hpp"

#include <stdexcept>

#include "vbscf/integrals/ao/one_electron/builder.hpp"

namespace xmvb::vb {

Eigen::MatrixXd ClosedShellFockBuilder::build(
    const Eigen::Ref<const Eigen::MatrixXd>& density_projector,
    const AoIntegralInput& ao_integral_input) const {
  // The closed-shell projector and inactive VB density have the same 2J-K
  // convention. Share the exact forward operator, including its MPI service.
  const int n_bf = ao_integral_input.n_basis_functions;
  if (n_bf <= 0 || density_projector.rows() != n_bf || density_projector.cols() != n_bf)
    throw std::invalid_argument("invalid closed-shell density dimensions");
  // The old eightfold traversal implicitly averaged both density directions.
  // Preserve that contract before the triangular-storage H1E operator.
  const Eigen::MatrixXd density = 0.5 * (density_projector + density_projector.transpose());
  return AoEffectiveOneElectronBuilder().build(
      density, ao_integral_input).ao_effective_h1e;
}

Eigen::MatrixXd ClosedShellFockBuilder::build(
    const Eigen::Ref<const Eigen::MatrixXd>& density_projector,
    const Eigen::Ref<const Eigen::MatrixXd>& core_hamiltonian_matrix,
    const RiAoFactorization& ri_factorization) const {
  const int n_basis_functions = ri_factorization.n_basis_functions;
  if (n_basis_functions <= 0 ||
      density_projector.rows() != n_basis_functions ||
      density_projector.cols() != n_basis_functions ||
      core_hamiltonian_matrix.rows() != n_basis_functions ||
      core_hamiltonian_matrix.cols() != n_basis_functions) {
    throw std::invalid_argument("invalid RI closed-shell Fock dimensions");
  }
  AoEffectiveOneElectronBuilder builder;
  return builder.build(
      density_projector,
      core_hamiltonian_matrix,
      ri_factorization,
      n_basis_functions).ao_effective_h1e;
}

}  // namespace xmvb::vb
