#include "integrals/active/ri.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <utility>
#include "core/pair_index.hpp"
#include "integrals/ao/ri_transform.hpp"
#include "parallel/openmp.hpp"

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

  ActiveSpaceTwoElectronResult result;
  result.representation =
      ActiveSpaceTwoElectronRepresentation::ResolutionOfIdentity;
  result.n_auxiliary_functions = ao_ri_result.n_auxiliary_functions;
  const int n_aux=ao_ri_result.n_auxiliary_functions;
  const int pairs=n_active_orbitals*(n_active_orbitals+1)/2;
  result.ri_active_pair_factors.resize(n_aux,pairs);
  const int workers=std::min(effective_openmp_thread_count(),n_basis_functions);
  const std::size_t per_aux=static_cast<std::size_t>(n_basis_functions)*
      (n_active_orbitals+workers)+n_active_orbitals*n_active_orbitals;
  const int tile=std::min<std::size_t>(n_aux,std::max<std::size_t>(1,
      ri_response_workspace_bytes/sizeof(double)/per_aux));
  std::shared_ptr<Eigen::MatrixXd> resident;
  if (static_cast<std::size_t>(n_aux)*n_active_orbitals*n_basis_functions<=
      ri_response_workspace_bytes/sizeof(double))
    resident=std::make_shared<Eigen::MatrixXd>(
        static_cast<Eigen::Index>(n_aux)*n_active_orbitals,n_basis_functions);
  for (int first=0;first<n_aux;first+=tile) {
    const int count=std::min(tile,n_aux-first);
    const Eigen::MatrixXd one_index=transform_ri_one_index(
        ao_ri_result.metric_whitened_ao_pair_factors,dense_active_coefficients,first,count);
    const Eigen::MatrixXd active=one_index*dense_active_coefficients;
    for (int a=0;a<n_active_orbitals;++a) {
      if (resident)
        resident->middleRows(a*n_aux+first,count)=one_index.middleRows(a*count,count);
      for (int b=0;b<=a;++b)
        result.ri_active_pair_factors.col(a*(a+1)/2+b).segment(first,count)=
            active.col(b).segment(a*count,count);
    }
  }
  result.ri_one_index_factors=std::move(resident);
  result.packed_active_two_electron_integrals =
      build_memory_dominated_pair_kernel(result.ri_active_pair_factors);
  result.dense_active_coefficients = dense_active_coefficients;

  return result;
}

}  // namespace xmvb::integrals
