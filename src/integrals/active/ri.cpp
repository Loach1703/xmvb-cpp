#include "integrals/active/ri.hpp"

#include <stdexcept>
#include "core/pair_index.hpp"
#include "integrals/active/pair_map.hpp"

namespace xmvb::integrals {
namespace {

std::vector<double> build_memory_dominated_pair_kernel(
    const Eigen::Ref<const Eigen::MatrixXd>& active_pair_factors) {
  const std::size_t n_pairs =
      static_cast<std::size_t>(active_pair_factors.cols());
  const std::size_t factor_elements =
      static_cast<std::size_t>(active_pair_factors.size());
  const std::size_t gram_elements = n_pairs * n_pairs;
  const std::size_t packed_elements = n_pairs * (n_pairs + 1) / 2;

  // A materialized pair kernel is useful only when its construction workspace
  // plus persistent packed result does not exceed one additional copy of the
  // already resident RI factor matrix. This dimension-only rule gives small
  // and medium active spaces one cache-friendly GEMM while preventing an
  // unconditional O(n_active^4) allocation for large active spaces.
  if (gram_elements + packed_elements > factor_elements) {
    return {};
  }

  const Eigen::MatrixXd gram =
      active_pair_factors.transpose() * active_pair_factors;
  std::vector<double> packed(packed_elements, 0.0);
  for (int column = 0; column < active_pair_factors.cols(); ++column) {
    for (int row = 0; row <= column; ++row) {
      packed[xmvb::core::PairIndex::packed_pair_of_pairs_index(row, column)] =
          gram(row, column);
    }
  }
  return packed;
}

}  // namespace

ActiveSpaceTwoElectronResult transform_ri_active_integrals(
    const RiAoFactorization& ao_ri_result,
    const Eigen::Ref<const Eigen::MatrixXd>& dense_active_coefficients) {
  const int n_basis_functions=static_cast<int>(dense_active_coefficients.rows());
  const int n_active_orbitals=static_cast<int>(dense_active_coefficients.cols());
  if (n_basis_functions <= 0 || n_active_orbitals <= 0) {
    throw std::invalid_argument("RI active-space dimensions must be positive");
  }
  if (ao_ri_result.n_basis_functions != n_basis_functions) {
    throw std::invalid_argument("AO RI result basis-function count mismatch");
  }

  const std::size_t n_ao_pairs =
      n_basis_functions * (n_basis_functions + 1) / 2;
  if (ao_ri_result.metric_whitened_ao_pair_factors.rows() !=
          ao_ri_result.n_auxiliary_functions ||
      ao_ri_result.metric_whitened_ao_pair_factors.cols() !=
          static_cast<Eigen::Index>(n_ao_pairs)) {
    throw std::invalid_argument("AO RI factor matrix shape mismatch");
  }

  PackedOrbitalPairMapMatrix ao_pair_to_active_pair_coefficients;
  build_packed_orbital_pair_map(
      dense_active_coefficients,
      &ao_pair_to_active_pair_coefficients);
  const Eigen::MatrixXd active_pair_factors =
      ao_ri_result.metric_whitened_ao_pair_factors *
      ao_pair_to_active_pair_coefficients;

  ActiveSpaceTwoElectronResult result;
  result.representation =
      ActiveSpaceTwoElectronRepresentation::ResolutionOfIdentity;
  result.n_auxiliary_functions = ao_ri_result.n_auxiliary_functions;
  result.ri_active_pair_factors = active_pair_factors;
  result.packed_active_two_electron_integrals =
      build_memory_dominated_pair_kernel(active_pair_factors);
  result.dense_active_coefficients = dense_active_coefficients;

  return result;
}

}  // namespace xmvb::integrals
