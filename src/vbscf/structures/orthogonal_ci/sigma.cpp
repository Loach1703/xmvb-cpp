#include "vbscf/structures/orthogonal_ci/sigma.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>

#include "core/openmp.hpp"
#include "vbscf/diagnostics/performance/phases.hpp"
#include "vbscf/integrals/active/two_electron/construction/indexer.hpp"
#include "vbscf/structures/orthogonal_ci/planner.hpp"

namespace xmvb::vb {
namespace {

std::uint64_t determinant_mask(
    const std::vector<int>& occupied,
    int n_orbitals) {
  std::uint64_t mask = 0;
  int previous = -1;
  for (const int orbital : occupied) {
    if (orbital <= previous || orbital < 0 || orbital >= n_orbitals) {
      throw std::invalid_argument("invalid direct-CI determinant occupation");
    }
    mask |= std::uint64_t{1} << orbital;
    previous = orbital;
  }
  return mask;
}

int only_set_bit(std::uint64_t mask) {
  if (mask == 0 || (mask & (mask - 1)) != 0) {
    throw std::logic_error("direct-CI single excitation is inconsistent");
  }
  return __builtin_ctzll(mask);
}

double replacement_sign(std::uint64_t source, int removed, int inserted) {
  const int low = std::min(removed, inserted);
  const int high = std::max(removed, inserted);
  std::uint64_t between = 0;
  if (high - low > 1) {
    between = source &
        (((std::uint64_t{1} << high) - 1) ^
         ((std::uint64_t{1} << (low + 1)) - 1));
  }
  return __builtin_popcountll(between) % 2 == 0 ? 1.0 : -1.0;
}

void add_symmetric_entry(
    Eigen::MatrixXd* matrix,
    int row,
    int column,
    double value) {
  if (row == column) {
    (*matrix)(row, column) += value;
    return;
  }
  (*matrix)(row, column) += 0.5 * value;
  (*matrix)(column, row) += 0.5 * value;
}

}  // namespace

void DirectCiSigmaAction::DensityConnections::reserve(std::size_t capacity) {
  sources.reserve(capacity);
  pairs.reserve(capacity);
  signs.reserve(capacity);
  created_orbitals.reserve(capacity);
  annihilated_orbitals.reserve(capacity);
}

void DirectCiSigmaAction::DensityConnections::append(
    int source,
    int pair,
    double sign,
    int created_orbital,
    int annihilated_orbital) {
  sources.push_back(source);
  pairs.push_back(pair);
  signs.push_back(sign);
  created_orbitals.push_back(created_orbital);
  annihilated_orbitals.push_back(annihilated_orbital);
}

std::size_t
DirectCiSigmaAction::DensityConnections::dynamic_bytes() const noexcept {
  return sources.capacity() * sizeof(int) +
      pairs.capacity() * sizeof(int) +
      signs.capacity() * sizeof(double) +
      created_orbitals.capacity() * sizeof(int) +
      annihilated_orbitals.capacity() * sizeof(int);
}

DirectCiSigmaAction::DirectCiSigmaAction(
    const std::vector<std::vector<int>>& alpha_determinants,
    const std::vector<std::vector<int>>& beta_determinants,
    const OrthogonalActiveIntegrals& integrals)
    : n_orbitals_(static_cast<int>(integrals.one_electron.rows())),
      n_alpha_(static_cast<int>(alpha_determinants.size())),
      n_beta_(static_cast<int>(beta_determinants.size())),
      pair_kernel_(integrals.pair_kernel) {
  if (n_orbitals_ <= 0 || n_orbitals_ > 63 ||
      integrals.one_electron.cols() != n_orbitals_ ||
      n_alpha_ <= 0 || n_beta_ <= 0) {
    throw std::invalid_argument("invalid orthogonal direct-CI sigma input");
  }
  if (!plan_orthogonal_direct_ci_action(
           alpha_determinants,
           beta_determinants,
           n_orbitals_,
           1).carrier_complete()) {
    throw std::invalid_argument(
        "orthogonal direct-CI sigma requires complete fixed-spin spaces");
  }
  const int n_electrons =
      static_cast<int>(alpha_determinants.front().size() +
                       beta_determinants.front().size());
  hamiltonian_pair_kernel_ =
      build_hamiltonian_pair_kernel(integrals, n_electrons);
  alpha_ = std::make_shared<const SpinConnections>(
      build_spin_connections(alpha_determinants));
  if (beta_determinants != alpha_determinants) {
    distinct_beta_ = std::make_shared<const SpinConnections>(
        build_spin_connections(beta_determinants));
  }
}

DirectCiSigmaAction::DirectCiSigmaAction(
    const DirectCiSigmaAction& reference,
    const OrthogonalActiveIntegrals& integrals)
    : n_orbitals_(reference.n_orbitals_),
      n_alpha_(reference.n_alpha_),
      n_beta_(reference.n_beta_),
      pair_kernel_(integrals.pair_kernel),
      alpha_(reference.alpha_),
      distinct_beta_(reference.distinct_beta_) {
  if (integrals.one_electron.rows() != n_orbitals_ ||
      integrals.pair_kernel.rows() != reference.pair_kernel_.rows()) {
    throw std::invalid_argument("direct-CI shared links have incompatible orbitals");
  }
  const int n_electrons = static_cast<int>(
      alpha_->occupied.front().size() + beta_connections().occupied.front().size());
  hamiltonian_pair_kernel_ = build_hamiltonian_pair_kernel(integrals, n_electrons);
}

DirectCiSigmaAction::SpinConnections
DirectCiSigmaAction::build_spin_connections(
    const std::vector<std::vector<int>>& determinants) const {
  diagnostics::ScopedPhase phase(diagnostics::Phase::CiConnections);
  const int dimension = static_cast<int>(determinants.size());
  std::vector<std::uint64_t> masks(dimension);
  std::unordered_map<std::uint64_t, int> index_by_mask;
  index_by_mask.reserve(determinants.size());
  for (int determinant = 0; determinant < dimension; ++determinant) {
    masks[determinant] = determinant_mask(
        determinants[determinant], n_orbitals_);
    if (!index_by_mask.emplace(masks[determinant], determinant).second) {
      throw std::invalid_argument(
          "direct-CI determinant space contains a duplicate");
    }
  }
  SpinConnections result;
  result.singles.resize(dimension);
  result.occupied = determinants;
  const int n_threads = coarse_openmp_thread_count(
      effective_openmp_thread_count(), static_cast<std::size_t>(dimension));

#pragma omp parallel for schedule(static) if(n_threads > 1) num_threads(n_threads)
  for (int target = 0; target < dimension; ++target) {
    auto& singles = result.singles[target];
    const auto& target_occupied = determinants[target];
    std::vector<int> target_virtual;
    target_virtual.reserve(n_orbitals_ - target_occupied.size());
    for (int orbital = 0; orbital < n_orbitals_; ++orbital) {
      if ((masks[target] & (std::uint64_t{1} << orbital)) == 0) {
        target_virtual.push_back(orbital);
      }
    }
    const std::size_t n_singles =
        target_occupied.size() * target_virtual.size();
    singles.reserve(n_singles);

    const auto add_source = [&](std::uint64_t source_mask) {
      const auto source_entry = index_by_mask.find(source_mask);
      if (source_entry == index_by_mask.end()) {
        throw std::logic_error(
            "direct-CI excitation left the fixed-spin determinant space");
      }
      const int source = source_entry->second;
      const int removed = only_set_bit(source_mask & ~masks[target]);
      const int inserted = only_set_bit(masks[target] & ~source_mask);
      const double sign = replacement_sign(source_mask, removed, inserted);
      singles.append(
          source,
          TwoElectronIndexer::packed_pair_index(removed, inserted),
          sign,
          inserted,
          removed);
    };

    for (const int removed_from_target : target_occupied) {
      for (const int inserted_into_source : target_virtual) {
        const std::uint64_t source_mask =
            (masks[target] ^
             (std::uint64_t{1} << removed_from_target)) |
            (std::uint64_t{1} << inserted_into_source);
        add_source(source_mask);
      }
    }
  }
  return result;
}

Eigen::MatrixXd DirectCiSigmaAction::build_hamiltonian_pair_kernel(
    const OrthogonalActiveIntegrals& integrals,
    int n_electrons) const {
  const int n_pairs = static_cast<int>(pair_kernel_.rows());
  if (n_electrons <= 0 || pair_kernel_.cols() != n_pairs ||
      integrals.one_electron.rows() != n_orbitals_ ||
      integrals.one_electron.cols() != n_orbitals_ ||
      !pair_kernel_.allFinite()) {
    throw std::invalid_argument(
        "invalid direct-CI packed Hamiltonian kernel input");
  }

  // In an orthonormal active basis,
  //   H = sum_pq h_pq E_pq
  //     + 1/2 sum_pqrs (pq|rs) (E_pq E_rs - delta_qr E_ps).
  // Absorb the one-body and contraction terms into the packed pair kernel so
  // the complete Hamiltonian is one scatter--GEMM--gather operation.
  Eigen::MatrixXd effective_one_electron = integrals.one_electron;
  for (int row = 0; row < n_orbitals_; ++row) {
    for (int column = 0; column < n_orbitals_; ++column) {
      double contraction = 0.0;
      for (int orbital = 0; orbital < n_orbitals_; ++orbital) {
        contraction += pair_kernel_(
            TwoElectronIndexer::packed_pair_index(row, orbital),
            TwoElectronIndexer::packed_pair_index(orbital, column));
      }
      effective_one_electron(row, column) -= 0.5 * contraction;
    }
  }
  effective_one_electron /= static_cast<double>(n_electrons);

  Eigen::MatrixXd result = 0.5 * pair_kernel_;
  for (int orbital = 0; orbital < n_orbitals_; ++orbital) {
    const int diagonal_pair =
        TwoElectronIndexer::packed_pair_index(orbital, orbital);
    for (int row = 0; row < n_orbitals_; ++row) {
      for (int column = 0; column <= row; ++column) {
        const int pair =
            TwoElectronIndexer::packed_pair_index(row, column);
        const double value = 0.5 * effective_one_electron(row, column);
        result(diagonal_pair, pair) += value;
        result(pair, diagonal_pair) += value;
      }
    }
  }
  return 0.5 * (result + result.transpose());
}

const DirectCiSigmaAction::SpinConnections&
DirectCiSigmaAction::beta_connections() const noexcept {
  return distinct_beta_ ? *distinct_beta_ : *alpha_;
}

void DirectCiSigmaAction::scatter_one_body_intermediate(
    const Eigen::Ref<const Eigen::MatrixXd>& coefficients,
    int alpha_begin,
    int alpha_count,
    int beta_target,
    Eigen::Ref<Eigen::MatrixXd> intermediate) const {
  const int block_width = static_cast<int>(coefficients.cols()) / n_beta_;
  const int active_columns = alpha_count * block_width;
  const SpinConnections& beta_graph = beta_connections();
  if (alpha_begin < 0 || alpha_count <= 0 ||
      alpha_begin + alpha_count > n_alpha_ ||
      beta_target < 0 || beta_target >= n_beta_ ||
      intermediate.rows() != active_columns ||
      intermediate.cols() != hamiltonian_pair_kernel_.rows()) {
    throw std::invalid_argument(
        "direct-CI scatter intermediate has incompatible dimensions");
  }
  intermediate.setZero();

  const auto add_beta_link = [&] (
      int beta_source, int pair, double sign) {
    for (int block = 0; block < block_width; ++block) {
      intermediate.col(pair)
          .segment(block * alpha_count, alpha_count)
          .noalias() += sign * coefficients.block(
              alpha_begin,
              block * n_beta_ + beta_source,
              alpha_count,
              1);
    }
  };
  for (const int occupied : beta_graph.occupied[beta_target]) {
    add_beta_link(
        beta_target,
        TwoElectronIndexer::packed_pair_index(occupied, occupied),
        1.0);
  }
  const DensityConnections& beta_singles = beta_graph.singles[beta_target];
  for (std::size_t link = 0; link < beta_singles.size(); ++link) {
    add_beta_link(
        beta_singles.sources[link],
        beta_singles.pairs[link],
        beta_singles.signs[link]);
  }

  for (int local_alpha = 0; local_alpha < alpha_count; ++local_alpha) {
    const int alpha_target = alpha_begin + local_alpha;
    const auto add_alpha_link = [&] (
        int alpha_source, int pair, double sign) {
      for (int block = 0; block < block_width; ++block) {
        intermediate(block * alpha_count + local_alpha, pair) +=
            sign * coefficients(
                alpha_source, block * n_beta_ + beta_target);
      }
    };
    for (const int occupied : alpha_->occupied[alpha_target]) {
      add_alpha_link(
          alpha_target,
          TwoElectronIndexer::packed_pair_index(occupied, occupied),
          1.0);
    }
    const DensityConnections& alpha_singles = alpha_->singles[alpha_target];
    for (std::size_t link = 0; link < alpha_singles.size(); ++link) {
      add_alpha_link(
          alpha_singles.sources[link],
          alpha_singles.pairs[link],
          alpha_singles.signs[link]);
    }
  }
}

Eigen::MatrixXd DirectCiSigmaAction::apply(
    const Eigen::Ref<const Eigen::MatrixXd>& coefficients) const {
  diagnostics::ScopedPhase phase(diagnostics::Phase::CiSigma);
  if (coefficients.rows() != n_alpha_ || coefficients.cols() <= 0 ||
      coefficients.cols() % n_beta_ != 0) {
    throw std::invalid_argument(
        "direct-CI coefficient block has incompatible dimensions");
  }
  const int block_width = static_cast<int>(coefficients.cols()) / n_beta_;
  const SpinConnections& beta_graph = beta_connections();
  const int n_pairs = static_cast<int>(hamiltonian_pair_kernel_.rows());
  Eigen::MatrixXd sigma = Eigen::MatrixXd::Zero(
      n_alpha_, coefficients.cols());
  const int n_threads = std::max(
      1, std::min(effective_openmp_thread_count(), n_beta_));

  // Bound all thread-local scatter, transformed-scatter, and reduction
  // buffers by one explicit workspace budget.  A wider alpha tile increases
  // GEMM efficiency without changing the mathematical action.
  constexpr std::size_t kWorkspaceBytes = 64ULL * 1024ULL * 1024ULL;
  constexpr int kMaximumAlphaTile = 160;
  const std::size_t bytes_per_alpha =
      static_cast<std::size_t>(n_threads) * block_width *
      (static_cast<std::size_t>(n_beta_) + 2ULL * n_pairs) * sizeof(double);
  const int alpha_tile = std::min(
      n_alpha_,
      std::max(
          1,
          std::min(
              kMaximumAlphaTile,
              static_cast<int>(kWorkspaceBytes /
                  std::max<std::size_t>(bytes_per_alpha, 1)))));
  const int maximum_columns = alpha_tile * block_width;
  std::vector<Eigen::MatrixXd> partial_sigma(n_threads);

#pragma omp parallel if(n_threads > 1) num_threads(n_threads)
  {
    int thread = 0;
    int team_size = 1;
#ifdef _OPENMP
    thread = omp_get_thread_num();
    team_size = omp_get_num_threads();
#endif
    // Strings are the contiguous axis, matching coefficient columns. Both
    // scatter/gather and BLAS now consume contiguous string segments.
    Eigen::MatrixXd intermediate(maximum_columns, n_pairs);
    Eigen::MatrixXd transformed(maximum_columns, n_pairs);
    Eigen::MatrixXd& beta_reduction = partial_sigma[thread];
    beta_reduction.resize(maximum_columns, n_beta_);

    for (int alpha_begin = 0;
         alpha_begin < n_alpha_;
         alpha_begin += alpha_tile) {
      const int alpha_count = std::min(alpha_tile, n_alpha_ - alpha_begin);
      const int active_columns = alpha_count * block_width;
      beta_reduction.topRows(active_columns).setZero();

#pragma omp for schedule(static)
      for (int beta_target = 0; beta_target < n_beta_; ++beta_target) {
        auto active_intermediate = intermediate.topRows(active_columns);
        scatter_one_body_intermediate(
            coefficients,
            alpha_begin,
            alpha_count,
            beta_target,
            active_intermediate);

        auto active_transformed = transformed.topRows(active_columns);
        active_transformed.noalias() =
            active_intermediate * hamiltonian_pair_kernel_.transpose();

        // Alpha links write one beta column owned by this loop iteration, so
        // no synchronization is required.
        for (int local_alpha = 0;
             local_alpha < alpha_count;
             ++local_alpha) {
          const int alpha_target = alpha_begin + local_alpha;
          const auto gather_alpha_link = [&] (
              int alpha_source, int pair, double sign) {
            for (int block = 0; block < block_width; ++block) {
              sigma(alpha_source, block * n_beta_ + beta_target) +=
                  sign * active_transformed(
                      block * alpha_count + local_alpha, pair);
            }
          };
          for (const int occupied : alpha_->occupied[alpha_target]) {
            gather_alpha_link(
                alpha_target,
                TwoElectronIndexer::packed_pair_index(occupied, occupied),
                1.0);
          }
        const DensityConnections& alpha_singles =
              alpha_->singles[alpha_target];
          for (std::size_t link = 0; link < alpha_singles.size(); ++link) {
            gather_alpha_link(
                alpha_singles.sources[link],
                alpha_singles.pairs[link],
                alpha_singles.signs[link]);
          }
        }

        // Contributions scattered through beta links can collide between
        // beta targets.  Accumulate them in one bounded buffer per thread.
        const auto gather_beta_link = [&] (
            int beta_source, int pair, double sign) {
          const double* source = active_transformed.col(pair).data();
          double* target = beta_reduction.col(beta_source).data();
#pragma omp simd
          for (int column = 0; column < active_columns; ++column) {
            target[column] += sign * source[column];
          }
        };
        for (const int occupied : beta_graph.occupied[beta_target]) {
          gather_beta_link(
              beta_target,
              TwoElectronIndexer::packed_pair_index(occupied, occupied),
              1.0);
        }
        const DensityConnections& beta_singles =
            beta_graph.singles[beta_target];
        for (std::size_t link = 0; link < beta_singles.size(); ++link) {
          gather_beta_link(
              beta_singles.sources[link],
              beta_singles.pairs[link],
              beta_singles.signs[link]);
        }
      }

#pragma omp for schedule(static)
      for (int beta = 0; beta < n_beta_; ++beta) {
        for (int block = 0; block < block_width; ++block) {
          for (int local_alpha = 0;
               local_alpha < alpha_count;
               ++local_alpha) {
            const int row = block * alpha_count + local_alpha;
            double value = 0.0;
            for (int source_thread = 0;
                 source_thread < team_size;
                 ++source_thread) {
              value += partial_sigma[source_thread](row, beta);
            }
            sigma(
                alpha_begin + local_alpha,
                block * n_beta_ + beta) += value;
          }
        }
      }
    }
  }
  return sigma;
}

Eigen::MatrixXd DirectCiSigmaAction::apply_one_body_generator(
    const Eigen::Ref<const Eigen::MatrixXd>& coefficients,
    const Eigen::Ref<const Eigen::MatrixXd>& generator) const {
  Eigen::MatrixXd result;
  apply_one_body_generators({{coefficients, generator, &result}});
  return result;
}

void DirectCiSigmaAction::apply_one_body_generators(
    std::initializer_list<DirectCiGeneratorRequest> requests) const {
  diagnostics::ScopedPhase phase(diagnostics::Phase::CiGenerator);
  if (requests.size() == 0) return;
  const Eigen::Index columns = requests.begin()->coefficients.cols();
  for (const auto& request : requests) {
    if (request.result == nullptr || request.coefficients.rows() != n_alpha_ ||
        columns <= 0 || columns % n_beta_ != 0 ||
        request.coefficients.cols() != columns ||
        request.generator.rows() != n_orbitals_ ||
        request.generator.cols() != n_orbitals_ ||
        !request.coefficients.allFinite() || !request.generator.allFinite()) {
      throw std::invalid_argument(
          "direct-CI generator requests have incompatible dimensions");
    }
    request.result->resize(n_alpha_, columns);
  }
  const SpinConnections& beta_graph = beta_connections();
  const Eigen::Index work_items = columns * static_cast<Eigen::Index>(n_alpha_);
  const int n_threads = coarse_openmp_thread_count(
      effective_openmp_thread_count(), static_cast<std::size_t>(work_items));

#pragma omp parallel for schedule(static) if(n_threads > 1) num_threads(n_threads)
  for (Eigen::Index work = 0; work < work_items; ++work) {
    const int alpha = work % n_alpha_;
    const Eigen::Index packed_column = work / n_alpha_;
    const int beta = packed_column % n_beta_;
    const Eigen::Index block = packed_column / n_beta_;
    const Eigen::Index column = packed_column;
    for (const auto& request : requests) {
      const auto& coefficients = request.coefficients;
      const auto& generator = request.generator;
      double value = 0.0;
      for (const int orbital : alpha_->occupied[alpha]) {
        value += generator(orbital, orbital) * coefficients(alpha, column);
      }
      for (const int orbital : beta_graph.occupied[beta]) {
        value += generator(orbital, orbital) * coefficients(alpha, column);
      }
      const DensityConnections& alpha_singles = alpha_->singles[alpha];
      for (std::size_t single = 0; single < alpha_singles.size(); ++single) {
        value += alpha_singles.signs[single] *
            generator(
                alpha_singles.created_orbitals[single],
                alpha_singles.annihilated_orbitals[single]) *
            coefficients(alpha_singles.sources[single], column);
      }
      const DensityConnections& beta_singles = beta_graph.singles[beta];
      for (std::size_t single = 0; single < beta_singles.size(); ++single) {
        value += beta_singles.signs[single] *
            generator(
                beta_singles.created_orbitals[single],
                beta_singles.annihilated_orbitals[single]) *
            coefficients(
                alpha,
                block * n_beta_ + beta_singles.sources[single]);
      }
      (*request.result)(alpha, column) = value;
    }
  }
}

DirectCiIntegralAdjoint DirectCiSigmaAction::integral_adjoint(
    const Eigen::Ref<const Eigen::MatrixXd>& left,
    const Eigen::Ref<const Eigen::MatrixXd>& right) const {
  if (left.rows() != n_alpha_ || right.rows() != n_alpha_ ||
      left.cols() <= 0 || left.cols() != right.cols() ||
      left.cols() % n_beta_ != 0 || !left.allFinite() ||
      !right.allFinite()) {
    throw std::invalid_argument(
        "direct-CI adjoint coefficient blocks have incompatible dimensions");
  }

  const int n_pairs = static_cast<int>(pair_kernel_.rows());
  const int block_width = static_cast<int>(left.cols()) / n_beta_;
  // Flatten the contracted spin-string rows into two shared, bounded tiles.
  // Each output block has one owner: no thread-private pair-square gradient
  // and no serial reduction of those gradients are needed.
  constexpr std::size_t workspace_bytes = 64ULL * 1024ULL * 1024ULL;
  const std::size_t bytes_per_row = 2ULL * n_pairs * sizeof(double);
  if (bytes_per_row > workspace_bytes) {
    throw std::runtime_error("one direct-CI adjoint row exceeds workspace budget");
  }
  const int block_tile = static_cast<int>(std::min<std::size_t>(
      block_width, workspace_bytes / bytes_per_row));
  const std::size_t bytes_per_alpha =
      bytes_per_row * block_tile;
  const int alpha_tile = static_cast<int>(std::min<std::size_t>(
      n_alpha_, workspace_bytes / bytes_per_alpha));
  const int rows_per_beta = alpha_tile * block_tile;
  const int beta_tile = static_cast<int>(std::min<std::size_t>(
      n_beta_, workspace_bytes / (bytes_per_alpha * alpha_tile)));
  const int maximum_rows = rows_per_beta * beta_tile;
  const int n_threads = std::max(
      1,
      std::min(effective_openmp_thread_count(), std::max(beta_tile, n_pairs)));
  // Among exact factorizations of the worker count, the most nearly square
  // output grid minimizes repeated reads of these equally sized operands.
  int row_parts = static_cast<int>(std::sqrt(n_threads));
  while (n_threads % row_parts != 0) --row_parts;
  const int column_parts = n_threads / row_parts;
  Eigen::MatrixXd left_intermediate(maximum_rows, n_pairs);
  Eigen::MatrixXd right_intermediate(maximum_rows, n_pairs);
  Eigen::MatrixXd effective_kernel_adjoint = Eigen::MatrixXd::Zero(n_pairs, n_pairs);

#pragma omp parallel if(n_threads > 1) num_threads(n_threads)
  {
    for (int block_begin = 0; block_begin < block_width; block_begin += block_tile) {
      const int block_count = std::min(block_tile, block_width - block_begin);
      const auto left_block = left.middleCols(block_begin * n_beta_, block_count * n_beta_);
      const auto right_block = right.middleCols(block_begin * n_beta_, block_count * n_beta_);
      for (int alpha_begin = 0; alpha_begin < n_alpha_; alpha_begin += alpha_tile) {
        const int alpha_count = std::min(alpha_tile, n_alpha_ - alpha_begin);
        const int rows = alpha_count * block_count;
        for (int beta_begin = 0; beta_begin < n_beta_; beta_begin += beta_tile) {
          const int beta_count = std::min(beta_tile, n_beta_ - beta_begin);
#pragma omp for schedule(static)
          for (int beta = 0; beta < beta_count; ++beta) {
            auto active_left = left_intermediate.middleRows(beta * rows, rows);
            auto active_right = right_intermediate.middleRows(beta * rows, rows);
            scatter_one_body_intermediate(
                left_block, alpha_begin, alpha_count, beta_begin + beta, active_left);
            scatter_one_body_intermediate(
                right_block, alpha_begin, alpha_count, beta_begin + beta, active_right);
          }
          // The scatter barrier publishes every shared row before GEMM; the
          // following barrier prevents overwriting the tile during contraction.
#pragma omp for schedule(static)
          for (int partition = 0; partition < n_threads; ++partition) {
            const int row_partition = partition / column_parts;
            const int column_partition = partition % column_parts;
            const int row_first = static_cast<int>(
                static_cast<std::int64_t>(n_pairs) * row_partition / row_parts);
            const int row_end = static_cast<int>(
                static_cast<std::int64_t>(n_pairs) * (row_partition + 1) / row_parts);
            const int first = static_cast<int>(
                static_cast<std::int64_t>(n_pairs) * column_partition / column_parts);
            const int end = static_cast<int>(
                static_cast<std::int64_t>(n_pairs) * (column_partition + 1) / column_parts);
            if (first == end || row_first == row_end) continue;
            effective_kernel_adjoint.block(
                row_first, first, row_end - row_first, end - first).noalias() +=
                left_intermediate.topRows(rows * beta_count)
                    .middleCols(row_first, row_end - row_first).transpose() *
                right_intermediate.topRows(rows * beta_count).middleCols(first, end - first);
          }
        }
      }
    }
  }

  Eigen::VectorXd absorbed_one_adjoint = Eigen::VectorXd::Zero(n_pairs);
  for (int orbital = 0; orbital < n_orbitals_; ++orbital) {
    const int diagonal_pair =
        TwoElectronIndexer::packed_pair_index(orbital, orbital);
    absorbed_one_adjoint.noalias() +=
        0.5 * (effective_kernel_adjoint.row(diagonal_pair).transpose() +
               effective_kernel_adjoint.col(diagonal_pair));
  }

  DirectCiIntegralAdjoint result;
  result.one_electron = Eigen::MatrixXd::Zero(
      n_orbitals_, n_orbitals_);
  result.pair_kernel =
      0.25 * (effective_kernel_adjoint +
              effective_kernel_adjoint.transpose());
  const int n_electrons =
      static_cast<int>(alpha_->occupied.front().size() +
                       beta_connections().occupied.front().size());
  for (int row = 0; row < n_orbitals_; ++row) {
    for (int column = 0; column <= row; ++column) {
      const int pair =
          TwoElectronIndexer::packed_pair_index(row, column);
      const double one_weight =
          absorbed_one_adjoint[pair] / static_cast<double>(n_electrons);
      add_symmetric_entry(
          &result.one_electron, row, column, one_weight);
      for (int orbital = 0; orbital < n_orbitals_; ++orbital) {
        add_symmetric_entry(
            &result.pair_kernel,
            TwoElectronIndexer::packed_pair_index(row, orbital),
            TwoElectronIndexer::packed_pair_index(orbital, column),
            -0.5 * one_weight);
      }
    }
  }
  return result;
}

Eigen::MatrixXd DirectCiSigmaAction::one_body_generator_adjoint(
    const Eigen::Ref<const Eigen::MatrixXd>& left,
    const Eigen::Ref<const Eigen::MatrixXd>& right) const {
  if (left.rows() != n_alpha_ || right.rows() != n_alpha_ ||
      left.cols() <= 0 || left.cols() != right.cols() ||
      left.cols() % n_beta_ != 0 || !left.allFinite() ||
      !right.allFinite()) {
    throw std::invalid_argument(
        "direct-CI generator-adjoint blocks have incompatible dimensions");
  }

  const int block_width = static_cast<int>(left.cols()) / n_beta_;
  const SpinConnections& beta_graph = beta_connections();
  const int work_items = block_width * n_alpha_ * n_beta_;
  const int n_threads = coarse_openmp_thread_count(
      effective_openmp_thread_count(), static_cast<std::size_t>(work_items));
  std::vector<Eigen::MatrixXd> partials(
      n_threads,
      Eigen::MatrixXd::Zero(n_orbitals_, n_orbitals_));

#pragma omp parallel if(n_threads > 1) num_threads(n_threads)
  {
    int thread = 0;
#ifdef _OPENMP
    thread = omp_get_thread_num();
#endif
    Eigen::MatrixXd& density = partials[thread];

#pragma omp for schedule(static)
    for (int work = 0; work < work_items; ++work) {
      const int beta = work % n_beta_;
      const int alpha = (work / n_beta_) % n_alpha_;
      const int block = work / (n_alpha_ * n_beta_);
      const int column = block * n_beta_ + beta;
      const double left_value = left(alpha, column);
      if (left_value == 0.0) {
        continue;
      }
      const double diagonal_weight =
          left_value * right(alpha, column);
      for (const int orbital : alpha_->occupied[alpha]) {
        density(orbital, orbital) += diagonal_weight;
      }
      for (const int orbital : beta_graph.occupied[beta]) {
        density(orbital, orbital) += diagonal_weight;
      }
      const DensityConnections& alpha_singles = alpha_->singles[alpha];
      for (std::size_t link = 0; link < alpha_singles.size(); ++link) {
        density(
            alpha_singles.annihilated_orbitals[link],
            alpha_singles.created_orbitals[link]) +=
            left_value * right(alpha_singles.sources[link], column) *
            alpha_singles.signs[link];
      }
      const DensityConnections& beta_singles = beta_graph.singles[beta];
      for (std::size_t link = 0; link < beta_singles.size(); ++link) {
        density(
            beta_singles.annihilated_orbitals[link],
            beta_singles.created_orbitals[link]) +=
            left_value * right(
                alpha,
                block * n_beta_ + beta_singles.sources[link]) *
            beta_singles.signs[link];
      }
    }
  }

  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(
      n_orbitals_, n_orbitals_);
  for (const Eigen::MatrixXd& partial : partials) {
    result += partial;
  }
  // The determinant Hamiltonian stores a one-body element for the
  // source-orbital/target-orbital pair. An orbital generator instead carries
  // the target/source convention in `E(kappa)`, hence the transpose.
  return result.transpose();
}

std::size_t DirectCiSigmaAction::dynamic_bytes() const noexcept {
  const auto spin_bytes = [](const SpinConnections& spin) {
    std::size_t bytes =
        spin.singles.capacity() * sizeof(DensityConnections) +
        spin.occupied.capacity() * sizeof(std::vector<int>);
    for (const DensityConnections& connections : spin.singles) {
      bytes += connections.dynamic_bytes();
    }
    for (const auto& occupied : spin.occupied) {
      bytes += occupied.capacity() * sizeof(int);
    }
    return bytes;
  };
  std::size_t bytes = static_cast<std::size_t>(
      pair_kernel_.size() + hamiltonian_pair_kernel_.size()) * sizeof(double) +
      spin_bytes(*alpha_);
  if (distinct_beta_) {
    bytes += spin_bytes(*distinct_beta_);
  }
  return bytes;
}

}  // namespace xmvb::vb
