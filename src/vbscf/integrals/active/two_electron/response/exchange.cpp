#include "vbscf/integrals/active/two_electron/response/internal.hpp"

#include <algorithm>
#include <stdexcept>

#include "core/openmp.hpp"
#include "vbscf/diagnostics/performance/phases.hpp"
#include "vbscf/integrals/ao/contracts/input.hpp"

namespace xmvb::vb::detail {

void build_exchange_integrals(ExactPackedActiveTwoElectronAdjointCache* cache) {
  const Eigen::Index n_bf = cache->n_basis_functions;
  const Eigen::Index n = cache->n_active_orbitals;
  constexpr std::size_t budget = 64 * 1024 * 1024;
  const std::size_t width = static_cast<std::size_t>(n) * n;
  const std::size_t rows = static_cast<std::size_t>(n_bf) * n_bf;
  const auto& graph = cache->ao_integral_input->pair_graph;
  // Keep the large-space streamed representation. The additional integral
  // image has its own memory bound and needs resident accepted pair products.
  if (cache->accepted_pair_products == nullptr || !graph.has_row_graph() ||
      width > budget / sizeof(double) / rows) return;
  // Include the shared active weights and mixed-direction workspace in the
  // admission bound, rather than bounding only the resident integral tensor.
  const std::size_t workspace = std::max(
      static_cast<std::size_t>(n_bf + n) * n * width,
      rows * static_cast<std::size_t>(n));
  if (rows * width + workspace > budget / sizeof(double)) return;
  diagnostics::ScopedPhase phase(diagnostics::Phase::ExchangeBuild);
  const ExactCtxPairMatrix coefficients = *cache->accepted_active_coefficients;
  cache->exchange_integrals.resize(n_bf * n_bf, n * n);
  const int threads = std::min<int>(effective_openmp_thread_count(), n_bf);
  const std::size_t slice_values = rows * static_cast<std::size_t>(n);
  const std::size_t available = budget / sizeof(double) - rows * width;
  const std::size_t tile_budget = std::min(available,
      std::max<std::size_t>(8 * 1024 * 1024 / sizeof(double), slice_values));
  const Eigen::Index tile = std::min<Eigen::Index>(n_bf, tile_budget / slice_values);
  ExactCtxPairMatrix half;
  // Transform one active index in the integral sweep, then the second with
  // GEMM. The half-transformed AO slices are shared and released per tile;
  // no worker owns a full image and the ERI loop has width n, not n*n.
  for (Eigen::Index begin = 0; begin < n_bf; begin += tile) {
    const Eigen::Index count = std::min(tile, n_bf - begin);
    half.resize(count * n_bf, n_bf * n);
#pragma omp parallel num_threads(threads)
    {
#pragma omp for schedule(dynamic, 1)
      for (Eigen::Index item = 0; item < count * n_bf; ++item) {
        const Eigen::Index mu = begin + item / n_bf;
        const Eigen::Index nu = item % n_bf;
        const auto pair = ao_pair_index(mu, nu);
        double* target = half.row(item).data();
        std::fill_n(target, n_bf * n, 0.0);
        for (int edge = graph.row_offsets[pair]; edge < graph.row_offsets[pair + 1]; ++edge) {
          const int right = graph.columns[edge];
          const int k = graph.pair_first[right], l = graph.pair_second[right];
          const double value = graph.values[edge];
          const double* c_k = coefficients.row(k).data();
          const double* c_l = coefficients.row(l).data();
#pragma omp simd
          for (Eigen::Index b = 0; b < n; ++b)
            target[l * n + b] += value * c_k[b];
          if (k != l) {
#pragma omp simd
            for (Eigen::Index b = 0; b < n; ++b)
              target[k * n + b] += value * c_l[b];
          }
        }
      }
#pragma omp for schedule(static)
      for (Eigen::Index item = 0; item < count * n_bf; ++item) {
        const Eigen::Index local_mu = item / n_bf;
        const Eigen::Index k = item % n_bf;
        Eigen::Map<ExactCtxPairMatrix> target(
            cache->exchange_integrals.row((begin + local_mu) * n_bf + k).data(), n, n);
        target.noalias() = coefficients.transpose() *
            half.middleRows(local_mu * n_bf, n_bf).middleCols(k * n, n);
      }
    }
  }
}

Eigen::MatrixXd contract_exchange_direction(
    const ExactPackedActiveTwoElectronAdjointCache& cache,
    const Eigen::Ref<const Eigen::MatrixXd>& direction) {
  diagnostics::ScopedPhase phase(diagnostics::Phase::ExchangeAction);
  const Eigen::Index n_bf = cache.n_basis_functions;
  const Eigen::Index n = cache.n_active_orbitals;
  if (direction.rows() != n_bf || direction.cols() != n ||
      cache.exchange_integrals.rows() != n_bf * n_bf ||
      cache.exchange_integrals.cols() != n * n)
    throw std::invalid_argument("exchange integral direction dimensions mismatch");

  // Expand packed pair multiplicities into a small active-only contraction.
  // Index (input orbital, partner-left, partner-right, output orbital).
  ExactCtxPairMatrix weights = ExactCtxPairMatrix::Zero(n, n * n * n);
  const Eigen::Index pairs = cache.active_pair_first_indices.size();
  for (Eigen::Index p = 0; p < pairs; ++p) {
    const int a = cache.active_pair_first_indices[p];
    const int b = cache.active_pair_second_indices[p];
    for (Eigen::Index q = 0; q < pairs; ++q) {
      const int c = cache.active_pair_first_indices[q];
      const int d = cache.active_pair_second_indices[q];
      const double value = cache.active_pair_gradient_matrix(q, p);
      weights(c, (b * n + d) * n + a) += value;
      weights(d, (b * n + c) * n + a) += value;
      weights(c, (a * n + d) * n + b) += value;
      weights(d, (a * n + c) * n + b) += value;
    }
  }
  const ExactCtxPairMatrix mixed = direction * weights;
  const Eigen::Map<const ExactCtxPairMatrix> source(mixed.data(), n_bf * n * n, n);
  const Eigen::Map<const ExactCtxPairMatrix> integrals(
      cache.exchange_integrals.data(), n_bf, n_bf * n * n);
  ExactCtxPairMatrix result(n_bf, n);
  const int threads = std::min<int>(effective_openmp_thread_count(), n_bf);
  const Eigen::Index tile = (n_bf + threads - 1) / threads;
#pragma omp parallel for num_threads(threads) schedule(static)
  for (Eigen::Index begin = 0; begin < n_bf; begin += tile) {
    const Eigen::Index count = std::min(tile, n_bf - begin);
    result.middleRows(begin, count).noalias() = integrals.middleRows(begin, count) * source;
  }
  return result;
}

}  // namespace xmvb::vb::detail
