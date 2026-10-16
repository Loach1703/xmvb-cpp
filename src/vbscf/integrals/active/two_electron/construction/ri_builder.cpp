#include "vbscf/integrals/active/two_electron/construction/ri_builder.hpp"

#include <stdexcept>
#include <vector>

#include <Eigen/Core>

#include "integrals/active/ri.hpp"

namespace xmvb::vb {

namespace {

Eigen::MatrixXd build_dense_active_coefficients(
    const OrbitalPreparationResult& orbital_preparation_result,
    int n_basis_functions,
    int n_active_orbitals) {
  if (orbital_preparation_result.active_sparse_row_offsets.size() !=
      n_basis_functions + 1) {
    throw std::invalid_argument("active_sparse_row_offsets size mismatch");
  }
  if (orbital_preparation_result.active_sparse_orbital_indices.size() !=
      orbital_preparation_result.active_sparse_values.size()) {
    throw std::invalid_argument("active sparse index/value sizes are inconsistent");
  }

  Eigen::MatrixXd dense_active_coefficients =
      Eigen::MatrixXd::Zero(n_basis_functions, n_active_orbitals);
  for (int basis_function_index = 0;
       basis_function_index < n_basis_functions;
       ++basis_function_index) {
    const int begin =
        orbital_preparation_result.active_sparse_row_offsets[
            basis_function_index];
    const int end =
        orbital_preparation_result.active_sparse_row_offsets[
            basis_function_index + 1];
    for (int offset = begin; offset < end; ++offset) {
      const int active_orbital_index =
          orbital_preparation_result.active_sparse_orbital_indices[
              offset];
      if (active_orbital_index < 0 || active_orbital_index >= n_active_orbitals) {
        throw std::invalid_argument("active sparse orbital index out of range");
      }
      dense_active_coefficients(basis_function_index, active_orbital_index) =
          orbital_preparation_result.active_sparse_values[offset];
    }
  }
  return dense_active_coefficients;
}

}  // namespace

integrals::ActiveSpaceTwoElectronResult RiActiveSpaceTwoElectronBuilder::build(
    const integrals::RiAoFactorization& ao,
    const OrbitalPreparationResult& orbitals,
    int n_bf, int n_active) const {
  if (n_bf<=0 || n_active<=0 || ao.n_basis_functions!=n_bf)
    throw std::invalid_argument("RI active-space dimensions mismatch");
  const auto active=build_dense_active_coefficients(orbitals,n_bf,n_active);
  return integrals::transform_ri_active_integrals(ao,active);
}

}  // namespace xmvb::vb
