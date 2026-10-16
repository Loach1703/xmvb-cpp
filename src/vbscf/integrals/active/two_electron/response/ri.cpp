#include "vbscf/integrals/active/two_electron/response/ri.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

#include <Eigen/Core>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "parallel/openmp.hpp"
#include "core/pair_index.hpp"

namespace xmvb::vb {
namespace {

std::size_t packed_pair_count(int dimension) {
  const std::size_t size = static_cast<std::size_t>(dimension);
  return size * (size + 1) / 2;
}

void validate_cache(
    const RiActiveTwoElectronResponseCache& cache,
    const char* operation) {
  const int n_bf = cache.n_basis_functions;
  const int n_active = cache.n_active_orbitals;
  const int n_auxiliary = cache.n_auxiliary_functions;
  if (n_bf <= 0 || n_active <= 0 || n_auxiliary <= 0) {
    throw std::invalid_argument(
        std::string(operation) + " requires positive RI dimensions");
  }
  if (cache.auxiliary_tile_size <= 0 ||
      cache.auxiliary_tile_size > n_auxiliary) {
    throw std::invalid_argument(
        std::string(operation) + " has an invalid auxiliary tile size");
  }
  const std::size_t n_ao_pairs = packed_pair_count(n_bf);
  const std::size_t n_active_pairs = packed_pair_count(n_active);
  if (cache.active_pair_first_indices.size() != n_active_pairs ||
      cache.active_pair_second_indices.size() != n_active_pairs) {
    throw std::invalid_argument(
        std::string(operation) + " has inconsistent pair-index tables");
  }
  if (cache.metric_whitened_ao_pair_factors == nullptr ||
      cache.metric_whitened_ao_pair_factors->rows() != n_auxiliary ||
      cache.metric_whitened_ao_pair_factors->cols() !=
          static_cast<Eigen::Index>(n_ao_pairs)) {
    throw std::invalid_argument(
        std::string(operation) + " has inconsistent AO RI factors");
  }
  if (cache.accepted_active_coefficients == nullptr ||
      cache.accepted_active_coefficients->rows() != n_bf ||
      cache.accepted_active_coefficients->cols() != n_active) {
    throw std::invalid_argument(
        std::string(operation) + " has inconsistent accepted coefficients");
  }
  if (cache.accepted_active_pair_factors == nullptr ||
      cache.accepted_active_pair_factors->rows() != n_auxiliary ||
      cache.accepted_active_pair_factors->cols() !=
          static_cast<Eigen::Index>(n_active_pairs)) {
    throw std::invalid_argument(
        std::string(operation) + " has inconsistent active-pair factors");
  }
  if (cache.accepted_one_index_factors.size() != 0 &&
      (cache.accepted_one_index_factors.rows() !=
           static_cast<Eigen::Index>(n_auxiliary) * n_active ||
       cache.accepted_one_index_factors.cols() != n_bf)) {
    throw std::invalid_argument(
        std::string(operation) + " has inconsistent one-index factors");
  }
}

Eigen::MatrixXd build_active_pair_adjoint(
    const std::vector<double>& packed_gradient,
    const RiActiveTwoElectronResponseCache& cache) {
  const Eigen::Index n_active_pairs = static_cast<Eigen::Index>(
      cache.active_pair_first_indices.size());
  const std::size_t expected_size =
      static_cast<std::size_t>(n_active_pairs) *
      (static_cast<std::size_t>(n_active_pairs) + 1) / 2;
  if (packed_gradient.size() != expected_size) {
    throw std::invalid_argument(
        "packed RI active-2e adjoint has the wrong size");
  }

  Eigen::MatrixXd active_pair_adjoint(n_active_pairs, n_active_pairs);
  for (Eigen::Index row = 0; row < n_active_pairs; ++row) {
    for (Eigen::Index column = 0; column < n_active_pairs; ++column) {
      const int first = static_cast<int>(row > column ? row : column);
      const int second = static_cast<int>(row > column ? column : row);
      const int packed_index =
          xmvb::core::PairIndex::packed_pair_of_pairs_index(first, second);
      double value = packed_gradient[packed_index];
      if (row == column) {
        value *= 2.0;
      }
      active_pair_adjoint(row, column) = value;
    }
  }
  return active_pair_adjoint;
}

/** @brief Shared factor/workspace budget; independent of numerical accuracy. */
constexpr std::size_t response_workspace_bytes = 64ULL * 1024 * 1024;

/** @brief Auxiliary tile width accounting for shared transforms and AO slices. */
int planned_auxiliary_tile_width(const RiActiveTwoElectronResponseCache& cache) {
  const std::size_t n_bf = cache.n_basis_functions;
  const std::size_t n_active = cache.n_active_orbitals;
  const std::size_t workers = std::min(
      xmvb::effective_openmp_thread_count(), cache.n_basis_functions);
  // Reserve two one-index transforms and one AO slice per worker.
  const std::size_t elements_per_aux = n_bf * (2 * n_active + workers);
  return static_cast<int>(std::min<std::size_t>(
      cache.n_auxiliary_functions,
      std::max<std::size_t>(
          1, response_workspace_bytes / sizeof(double) / elements_per_aux)));
}

/** @brief Forms L_A C without materializing full symmetric AO factors.
 *
 * Rows are grouped by active orbital, with auxiliary indices contiguous.
 * Workers own only one tiled AO slice; the output is shared.
 */
Eigen::MatrixXd transform_one_index(
    const RiActiveTwoElectronResponseCache& cache,
    const Eigen::Ref<const Eigen::MatrixXd>& coefficients,
    int first_aux, int count) {
  const int n_bf = cache.n_basis_functions;
  const int n_active = cache.n_active_orbitals;
  Eigen::MatrixXd transformed(
      static_cast<Eigen::Index>(count) * n_active, n_bf);
  const int workers = std::min(
      xmvb::effective_openmp_thread_count(), n_bf);
#pragma omp parallel num_threads(workers)
  {
    Eigen::MatrixXd ao_slice(count, n_bf);
#pragma omp for schedule(static)
    for (int mu = 0; mu < n_bf; ++mu) {
      for (int nu = 0; nu < n_bf; ++nu) {
        const int larger = std::max(mu, nu);
        const int smaller = std::min(mu, nu);
        const Eigen::Index pair =
            static_cast<Eigen::Index>(larger) * (larger + 1) / 2 + smaller;
        ao_slice.col(nu) =
            cache.metric_whitened_ao_pair_factors->col(pair).segment(
                first_aux, count);
      }
      Eigen::Map<Eigen::MatrixXd> output(
          transformed.col(mu).data(), count, n_active);
      output.noalias() = ao_slice * coefficients;
    }
  }
  return transformed;
}

/** @brief Visits shared accepted factors, or generates the same tiled factors. */
template <class Consumer>
void visit_accepted_factors(
    const RiActiveTwoElectronResponseCache& cache, Consumer consume) {
  if (cache.accepted_one_index_factors.size() != 0) {
    // No copy and no direction-dependent accepted transform.
    consume(cache.accepted_one_index_factors, 0, cache.n_auxiliary_functions);
    return;
  }
  const int width = cache.auxiliary_tile_size;
  for (int first = 0; first < cache.n_auxiliary_functions; first += width) {
    const int count = std::min(width, cache.n_auxiliary_functions - first);
    const Eigen::MatrixXd factors = transform_one_index(
        cache, *cache.accepted_active_coefficients, first, count);
    consume(factors, first, count);
  }
}

/** @brief delta B_ab = (L C_a)^T D_b + (L C_b)^T D_a. */
void factor_direction_tile(
    const RiActiveTwoElectronResponseCache& cache,
    const Eigen::MatrixXd& factors,
    const Eigen::Ref<const Eigen::MatrixXd>& direction,
    int first_aux, int count, Eigen::MatrixXd& result) {
  const Eigen::MatrixXd mixed = factors * direction;
  for (Eigen::Index pair = 0; pair < result.cols(); ++pair) {
    const int a = cache.active_pair_first_indices[pair];
    const int b = cache.active_pair_second_indices[pair];
    result.col(pair).segment(first_aux, count) =
        mixed.col(b).segment(a * count, count) +
        mixed.col(a).segment(b * count, count);
  }
}

/** @brief Unpacks a factor adjoint with canonical diagonal multiplicity. */
void unpack_factor_adjoint(
    const RiActiveTwoElectronResponseCache& cache,
    const Eigen::MatrixXd& adjoint, int auxiliary, Eigen::MatrixXd& weights) {
  for (Eigen::Index pair = 0; pair < adjoint.cols(); ++pair) {
    const int a = cache.active_pair_first_indices[pair];
    const int b = cache.active_pair_second_indices[pair];
    const double value = adjoint(auxiliary, pair);
    weights(a, b) = a == b ? 2.0 * value : value;
    weights(b, a) = weights(a, b);
  }
}

/** @brief Contracts factor adjoints directly to AO-by-active orbital output.
 *
 * The optional directional term differentiates L_A C, not an AO-pair
 * adjoint. Per-worker storage is AO-by-active plus active-by-active only.
 */
Eigen::MatrixXd contract_factor_adjoint(
    const RiActiveTwoElectronResponseCache& cache,
    const Eigen::MatrixXd& adjoint,
    const Eigen::MatrixXd* directional_adjoint,
    const Eigen::Ref<const Eigen::MatrixXd>* direction) {
  const int n_bf = cache.n_basis_functions;
  const int n_active = cache.n_active_orbitals;
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n_bf, n_active);
  visit_accepted_factors(cache, [&](const Eigen::MatrixXd& factors,
                                  int first_aux, int count) {
    // Repack only the one-index directional factors, never AO-pair adjoints.
    Eigen::MatrixXd directional_factors;
    const int width = cache.auxiliary_tile_size;
    for (int offset = 0; offset < count; offset += width) {
      const int tile_count = std::min(width, count - offset);
      if (direction != nullptr) {
        directional_factors = transform_one_index(
            cache, *direction, first_aux + offset, tile_count);
      }
      const std::size_t worker_elements =
          3ULL * n_bf * n_active + 2ULL * n_active * n_active;
      const int workers = xmvb::bounded_openmp_thread_count(
          static_cast<int>(std::min<std::size_t>(
              xmvb::effective_openmp_thread_count(),
              std::max<std::size_t>(
                  1, response_workspace_bytes / sizeof(double) /
                         worker_elements))),
          tile_count);
      std::vector<Eigen::MatrixXd> partial(
          workers, Eigen::MatrixXd::Zero(n_bf, n_active));
#pragma omp parallel num_threads(workers)
      {
        int worker = 0;
#ifdef _OPENMP
        worker = omp_get_thread_num();
#endif
        Eigen::MatrixXd image(n_bf, n_active);
        Eigen::MatrixXd weights(n_active, n_active);
        Eigen::MatrixXd delta_image(n_bf, n_active);
        Eigen::MatrixXd delta_weights(n_active, n_active);
#pragma omp for schedule(static)
        for (int auxiliary = 0; auxiliary < tile_count; ++auxiliary) {
          for (int a = 0; a < n_active; ++a) {
            image.col(a) = factors.row(
                a * count + offset + auxiliary).transpose();
          }
          unpack_factor_adjoint(
              cache, adjoint, first_aux + offset + auxiliary, weights);
          if (direction == nullptr) {
            partial[worker].noalias() += image * weights;
          } else {
            for (int a = 0; a < n_active; ++a) {
              delta_image.col(a) = directional_factors.row(
                  a * tile_count + auxiliary).transpose();
            }
            unpack_factor_adjoint(cache, *directional_adjoint,
                                  first_aux + offset + auxiliary,
                                  delta_weights);
            partial[worker].noalias() += delta_image * weights;
            partial[worker].noalias() += image * delta_weights;
          }
        }
      }
      for (const auto& value : partial) result += value;
    }
  });
  return result;
}
}  // namespace

RiActiveTwoElectronResponseCache build_ri_active_two_electron_response_cache(
    const RiAoFactorization& ao_ri_factorization,
    const ActiveSpaceTwoElectronResult& accepted_active_space_result,
    int n_active_orbitals) {
  if (ao_ri_factorization.n_basis_functions <= 0 ||
      ao_ri_factorization.n_auxiliary_functions <= 0 ||
      n_active_orbitals <= 0) {
    throw std::invalid_argument(
        "RI active-2e response cache dimensions must be positive");
  }
  if (accepted_active_space_result.representation !=
      ActiveSpaceTwoElectronRepresentation::ResolutionOfIdentity) {
    throw std::invalid_argument(
        "RI active-2e response cache requires an RI active-space result");
  }

  const int n_bf = ao_ri_factorization.n_basis_functions;
  const std::size_t n_ao_pairs = packed_pair_count(n_bf);
  const std::size_t n_active_pairs = packed_pair_count(n_active_orbitals);
  if (ao_ri_factorization.n_packed_ao_pairs !=
          static_cast<int>(n_ao_pairs) ||
      ao_ri_factorization.metric_whitened_ao_pair_factors.rows() !=
          ao_ri_factorization.n_auxiliary_functions ||
      ao_ri_factorization.metric_whitened_ao_pair_factors.cols() !=
          static_cast<Eigen::Index>(n_ao_pairs)) {
    throw std::invalid_argument(
        "RI active-2e response cache AO factor dimensions are inconsistent");
  }
  if (accepted_active_space_result.n_auxiliary_functions !=
          ao_ri_factorization.n_auxiliary_functions ||
      accepted_active_space_result.dense_active_coefficients.rows() != n_bf ||
      accepted_active_space_result.dense_active_coefficients.cols() !=
          n_active_orbitals ||
      accepted_active_space_result.ri_active_pair_factors.rows() !=
          ao_ri_factorization.n_auxiliary_functions ||
      accepted_active_space_result.ri_active_pair_factors.cols() !=
          static_cast<Eigen::Index>(n_active_pairs)) {
    throw std::invalid_argument(
        "RI active-2e response cache accepted dimensions are inconsistent");
  }

  RiActiveTwoElectronResponseCache cache;
  cache.n_basis_functions = n_bf;
  cache.n_active_orbitals = n_active_orbitals;
  cache.n_auxiliary_functions =
      ao_ri_factorization.n_auxiliary_functions;
  cache.metric_whitened_ao_pair_factors =
      &ao_ri_factorization.metric_whitened_ao_pair_factors;
  cache.accepted_active_coefficients =
      &accepted_active_space_result.dense_active_coefficients;
  cache.accepted_active_pair_factors =
      &accepted_active_space_result.ri_active_pair_factors;

  cache.active_pair_first_indices.reserve(n_active_pairs);
  cache.active_pair_second_indices.reserve(n_active_pairs);
  for (int first = 0; first < n_active_orbitals; ++first) {
    for (int second = 0; second <= first; ++second) {
      cache.active_pair_first_indices.push_back(first);
      cache.active_pair_second_indices.push_back(second);
    }
  }
  const std::size_t factor_elements =
      static_cast<std::size_t>(cache.n_auxiliary_functions) * n_bf *
      n_active_orbitals;
  cache.auxiliary_tile_size = planned_auxiliary_tile_width(cache);
  if (factor_elements <= response_workspace_bytes / sizeof(double)) {
    cache.accepted_one_index_factors.resize(
        static_cast<Eigen::Index>(cache.n_auxiliary_functions) *
            n_active_orbitals,
        n_bf);
    const int width = cache.auxiliary_tile_size;
    for (int first = 0; first < cache.n_auxiliary_functions; first += width) {
      const int count = std::min(width, cache.n_auxiliary_functions - first);
      const Eigen::MatrixXd factors = transform_one_index(
          cache, *cache.accepted_active_coefficients, first, count);
      for (int a = 0; a < n_active_orbitals; ++a) {
        cache.accepted_one_index_factors.middleRows(
            static_cast<Eigen::Index>(a) * cache.n_auxiliary_functions + first,
            count) = factors.middleRows(a * count, count);
      }
    }
  }
  validate_cache(cache, "RI active-2e response cache");
  return cache;
}

Eigen::MatrixXd compute_ri_active_pair_factor_directional_derivative(
    const RiActiveTwoElectronResponseCache& accepted_cache,
    const Eigen::Ref<const Eigen::MatrixXd>& dense_active_direction) {
  validate_cache(accepted_cache, "RI active-pair factor direction");
  if (dense_active_direction.rows() != accepted_cache.n_basis_functions ||
      dense_active_direction.cols() != accepted_cache.n_active_orbitals) {
    throw std::invalid_argument(
        "RI active-pair factor direction has the wrong shape");
  }
  Eigen::MatrixXd result(
      accepted_cache.n_auxiliary_functions,
      accepted_cache.active_pair_first_indices.size());
  visit_accepted_factors(accepted_cache,
      [&](const Eigen::MatrixXd& factors, int first, int count) {
        factor_direction_tile(accepted_cache, factors, dense_active_direction,
                              first, count, result);
      });
  return result;
}

std::vector<Eigen::MatrixXd>
compute_ri_active_pair_factor_directional_derivative_batch(
    const RiActiveTwoElectronResponseCache& accepted_cache,
    const std::vector<Eigen::MatrixXd>& dense_active_directions) {
  validate_cache(accepted_cache, "RI active-pair factor direction block");
  if (dense_active_directions.empty()) return {};
  for (const auto& direction : dense_active_directions) {
    if (direction.rows() != accepted_cache.n_basis_functions ||
        direction.cols() != accepted_cache.n_active_orbitals ||
        !direction.allFinite()) {
      throw std::invalid_argument(
          "RI active-pair factor direction block has an invalid shape or value");
    }
  }
  std::vector<Eigen::MatrixXd> results;
  results.reserve(dense_active_directions.size());
  for (std::size_t i = 0; i < dense_active_directions.size(); ++i) {
    results.emplace_back(
        accepted_cache.n_auxiliary_functions,
        accepted_cache.active_pair_first_indices.size());
  }
  visit_accepted_factors(accepted_cache,
      [&](const Eigen::MatrixXd& factors, int first, int count) {
        for (std::size_t i = 0; i < dense_active_directions.size(); ++i) {
          factor_direction_tile(accepted_cache, factors,
                                dense_active_directions[i], first, count,
                                results[i]);
        }
      });
  return results;
}

std::vector<double>
compute_ri_packed_active_two_electron_integral_directional_derivative(
    const RiActiveTwoElectronResponseCache& accepted_cache,
    const Eigen::Ref<const Eigen::MatrixXd>&
        directional_active_pair_factors) {
  validate_cache(accepted_cache, "RI packed active-2e direction");
  const Eigen::MatrixXd& accepted_factors =
      *accepted_cache.accepted_active_pair_factors;
  if (directional_active_pair_factors.rows() != accepted_factors.rows() ||
      directional_active_pair_factors.cols() != accepted_factors.cols()) {
    throw std::invalid_argument(
        "RI active-pair factor direction has the wrong shape");
  }

  const Eigen::MatrixXd one_sided_contraction =
      accepted_factors.transpose() * directional_active_pair_factors;
  const Eigen::Index n_active_pairs = accepted_factors.cols();
  std::vector<double> packed_direction(
      static_cast<std::size_t>(n_active_pairs) *
          (static_cast<std::size_t>(n_active_pairs) + 1) / 2,
      0.0);
  for (Eigen::Index first = 0; first < n_active_pairs; ++first) {
    for (Eigen::Index second = 0; second <= first; ++second) {
      const int packed_index =
          xmvb::core::PairIndex::packed_pair_of_pairs_index(
              static_cast<int>(first),
              static_cast<int>(second));
      packed_direction[packed_index] =
          one_sided_contraction(first, second) +
          one_sided_contraction(second, first);
    }
  }
  return packed_direction;
}

Eigen::MatrixXd
compute_ri_packed_active_two_electron_integral_directional_derivative_batch(
    const RiActiveTwoElectronResponseCache& accepted_cache,
    const std::vector<Eigen::MatrixXd>& directional_active_pair_factors) {
  validate_cache(accepted_cache, "RI packed active-2e direction block");
  if (directional_active_pair_factors.empty()) return {};
  const Eigen::MatrixXd& accepted_factors =
      *accepted_cache.accepted_active_pair_factors;
  const Eigen::Index n_active_pairs = accepted_factors.cols();
  const Eigen::Index n_directions = static_cast<Eigen::Index>(
      directional_active_pair_factors.size());
  Eigen::MatrixXd factor_direction_block(
      accepted_factors.rows(), n_directions * n_active_pairs);
  for (Eigen::Index direction = 0;
       direction < n_directions; ++direction) {
    const auto& factor_direction =
        directional_active_pair_factors[static_cast<std::size_t>(direction)];
    if (factor_direction.rows() != accepted_factors.rows() ||
        factor_direction.cols() != n_active_pairs ||
        !factor_direction.allFinite()) {
      throw std::invalid_argument(
          "RI active-pair factor direction block has an invalid shape or value");
    }
    factor_direction_block.middleCols(
        direction * n_active_pairs, n_active_pairs) = factor_direction;
  }
  const Eigen::MatrixXd one_sided_block =
      accepted_factors.transpose() * factor_direction_block;
  const Eigen::Index n_packed_integrals =
      n_active_pairs * (n_active_pairs + 1) / 2;
  Eigen::MatrixXd packed_directions(n_packed_integrals, n_directions);
  for (Eigen::Index direction = 0;
       direction < n_directions; ++direction) {
    const auto one_sided = one_sided_block.middleCols(
        direction * n_active_pairs, n_active_pairs);
    for (Eigen::Index first = 0; first < n_active_pairs; ++first) {
      for (Eigen::Index second = 0; second <= first; ++second) {
        const int packed_index =
            xmvb::core::PairIndex::packed_pair_of_pairs_index(
                static_cast<int>(first), static_cast<int>(second));
        packed_directions(packed_index, direction) =
            one_sided(first, second) + one_sided(second, first);
      }
    }
  }
  return packed_directions;
}

Eigen::MatrixXd backpropagate_ri_packed_active_two_electron_gradient(
    const std::vector<double>& packed_active_two_electron_gradient,
    const RiActiveTwoElectronResponseCache& accepted_cache) {
  validate_cache(accepted_cache, "RI packed active-2e pullback");
  const Eigen::MatrixXd active_pair_adjoint =
      build_active_pair_adjoint(
          packed_active_two_electron_gradient,
          accepted_cache);
  const Eigen::MatrixXd factor_adjoint =
      *accepted_cache.accepted_active_pair_factors * active_pair_adjoint;
  return contract_factor_adjoint(
      accepted_cache, factor_adjoint, nullptr, nullptr);
}

Eigen::MatrixXd
apply_ri_packed_active_two_electron_adjoint_hessian_vector(
    const std::vector<double>& packed_active_two_electron_gradient,
    const RiActiveTwoElectronResponseCache& accepted_cache,
    const Eigen::Ref<const Eigen::MatrixXd>& dense_active_direction,
    const Eigen::MatrixXd* directional_active_pair_factors) {
  validate_cache(accepted_cache, "RI fixed-adjoint active-2e HVP");
  if (dense_active_direction.rows() != accepted_cache.n_basis_functions ||
      dense_active_direction.cols() != accepted_cache.n_active_orbitals) {
    throw std::invalid_argument(
        "RI fixed-adjoint HVP direction has the wrong orbital shape");
  }

  Eigen::MatrixXd owned_factor_direction;
  const Eigen::MatrixXd* factor_direction =
      directional_active_pair_factors;
  if (factor_direction == nullptr) {
    owned_factor_direction =
        compute_ri_active_pair_factor_directional_derivative(
            accepted_cache,
            dense_active_direction);
    factor_direction = &owned_factor_direction;
  }
  if (factor_direction->rows() != accepted_cache.n_auxiliary_functions ||
      factor_direction->cols() !=
          static_cast<Eigen::Index>(
              accepted_cache.active_pair_first_indices.size())) {
    throw std::invalid_argument(
        "precomputed RI active-pair factor direction has the wrong shape");
  }

  const Eigen::MatrixXd active_pair_adjoint =
      build_active_pair_adjoint(
          packed_active_two_electron_gradient,
          accepted_cache);
  const Eigen::MatrixXd accepted_factor_adjoint =
      *accepted_cache.accepted_active_pair_factors * active_pair_adjoint;
  const Eigen::MatrixXd directional_factor_adjoint =
      *factor_direction * active_pair_adjoint;
  return contract_factor_adjoint(
      accepted_cache, accepted_factor_adjoint, &directional_factor_adjoint,
      &dense_active_direction);
}

}  // namespace xmvb::vb
