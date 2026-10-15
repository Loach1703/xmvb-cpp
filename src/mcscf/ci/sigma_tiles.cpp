#include "mcscf/ci/sigma.hpp"

#include <algorithm>
#include <exception>
#include <limits>
#include <stdexcept>

#include "core/openmp.hpp"
#include "core/pair_index.hpp"
#include "core/timing.hpp"
#include "mcscf/ci/covalent_support.hpp"

namespace xmvb::mcscf {

void DirectCiSigmaAction::apply_tiled(
    int block_width,
    const std::function<void(int, const std::vector<int>&, Eigen::Ref<Eigen::MatrixXd>)>& read_column,
    const std::function<void(int, const std::vector<int>&, const Eigen::Ref<const Eigen::MatrixXd>&)>& add_column,
    const CovalentSupport* support) const {
  xmvb::profiling::ScopedPhase phase(xmvb::profiling::Phase::CiSigma);
  if (block_width <= 0 || !read_column || !add_column ||
      static_cast<Eigen::Index>(block_width) >
          std::numeric_limits<Eigen::Index>::max() / n_alpha_) {
    throw std::invalid_argument("invalid tiled direct-CI sigma callbacks or width");
  }

  constexpr int kMaximumAlphaTile = 160;
  constexpr std::size_t kWorkspaceBytes = 64ULL * 1024ULL * 1024ULL;
  const int n_pairs = static_cast<int>(hamiltonian_pair_kernel_.rows());
  const SpinConnections& beta_graph = beta_connections();
  const int available_threads =
      std::max(1, std::min(effective_openmp_thread_count(), n_beta_));
  const auto bounded_product = [](std::size_t left, std::size_t right) {
    return right != 0 && left > std::numeric_limits<std::size_t>::max() / right
        ? std::numeric_limits<std::size_t>::max() : left * right;
  };

  std::vector<int> halo;
  std::vector<int> halo_index(n_alpha_, -1);
  const auto collect_halo = [&](int alpha_begin, int alpha_count) {
    for (const int alpha : halo) halo_index[alpha] = -1;
    halo.clear();
    const auto include = [&](int alpha) {
      if (halo_index[alpha] < 0) {
        halo_index[alpha] = 0;
        halo.push_back(alpha);
      }
    };
    for (int alpha_target = alpha_begin;
         alpha_target < alpha_begin + alpha_count; ++alpha_target) {
      include(alpha_target);
      const DensityConnections& links = alpha_->singles[alpha_target];
      for (const int source : links.sources) include(source);
    }
    std::sort(halo.begin(), halo.end());
    for (int local = 0; local < static_cast<int>(halo.size()); ++local)
      halo_index[halo[local]] = local;
  };

  // Plan one fixed team and workspace shape. Inspect link metadata only;
  // never construct a carrier or reject a physical space based on this budget.
  int alpha_extent = std::min(kMaximumAlphaTile, n_alpha_);
  std::size_t maximum_halo = 0, shared_bytes = 0, worker_bytes = 0;
  for (;;) {
    maximum_halo = 0;
    for (int first = 0; first < n_alpha_; first += alpha_extent) {
      collect_halo(first, std::min(alpha_extent, n_alpha_ - first));
      maximum_halo = std::max(maximum_halo, halo.size());
    }
    shared_bytes = bounded_product(2 * sizeof(double) + (support ? 1 : 0), bounded_product(
        bounded_product(maximum_halo, static_cast<std::size_t>(block_width)),
        static_cast<std::size_t>(n_beta_) + 1));
    worker_bytes = bounded_product(
          bounded_product(static_cast<std::size_t>(alpha_extent),
                          static_cast<std::size_t>(block_width)),
          bounded_product(static_cast<std::size_t>(n_beta_) + 2ULL * n_pairs,
                          sizeof(double)));
    if (alpha_extent == 1 ||
        (shared_bytes < kWorkspaceBytes &&
         worker_bytes <= (kWorkspaceBytes - shared_bytes) /
             static_cast<std::size_t>(available_threads))) break;
    alpha_extent = std::max(1, alpha_extent / 2);
  }

  const int n_threads = static_cast<int>(std::min<std::size_t>(
        static_cast<std::size_t>(available_threads),
        std::max<std::size_t>(1,
            shared_bytes >= kWorkspaceBytes ? 1 :
                (kWorkspaceBytes - shared_bytes) / worker_bytes)));
  const Eigen::Index halo_columns =
        static_cast<Eigen::Index>(block_width) * n_beta_;
  const int maximum_rows = alpha_extent * block_width;
  Eigen::MatrixXd input(maximum_halo, halo_columns), output(maximum_halo, halo_columns);
  Eigen::MatrixXd source_column(maximum_halo, block_width), output_column(maximum_halo, block_width);
  // This is a bounded halo tile, not a retained U_alpha*U_beta mask.
  std::vector<std::uint8_t> allowed(support ? maximum_halo * n_beta_ : 0);
  std::vector<int> column_rows;
  column_rows.reserve(maximum_halo);
  std::vector<Eigen::MatrixXd> intermediates(n_threads), transformed_tiles(n_threads), reductions(n_threads);
  std::vector<double*> partial_sigma(n_threads);
  // Allocate once before entering OpenMP, so allocation failures propagate
  // normally. Workers take ownership and first-touch their private buffers.
  for (int worker = 0; worker < n_threads; ++worker) {
    intermediates[worker].resize(maximum_rows, n_pairs);
    transformed_tiles[worker].resize(maximum_rows, n_pairs);
    reductions[worker].resize(maximum_rows, n_beta_);
    partial_sigma[worker] = reductions[worker].data();
  }
  std::exception_ptr failure;
#pragma omp parallel if(n_threads > 1) num_threads(n_threads)
  {
    int thread = 0, team_size = 1;
#ifdef _OPENMP
    thread = omp_get_thread_num();
    team_size = omp_get_num_threads();
#endif
    Eigen::MatrixXd intermediate = std::move(intermediates[thread]);
    Eigen::MatrixXd transformed = std::move(transformed_tiles[thread]);
    Eigen::MatrixXd beta_reduction = std::move(reductions[thread]);
    std::vector<int> active_alpha;
    active_alpha.reserve(alpha_extent);
    for (int alpha_begin = 0; alpha_begin < n_alpha_; alpha_begin += alpha_extent) {
      const int alpha_count = std::min(alpha_extent, n_alpha_ - alpha_begin);
      const int active_rows = alpha_count * block_width;
#pragma omp single
      {
        try {
          collect_halo(alpha_begin, alpha_count);
          output.topRows(halo.size()).setZero();
        } catch (...) { failure = std::current_exception(); }
      }
      if (failure) break;
      // Prefix tests are tile metadata, shared by both gather and adjoint.
      // Do not serialize this exact integer classification in the callbacks.
#pragma omp for schedule(static)
      for (int beta = 0; beta < n_beta_; ++beta)
        if (support)
          for (int local = 0; local < static_cast<int>(halo.size()); ++local)
            allowed[static_cast<std::size_t>(beta) * maximum_halo + local] =
                support->contains(halo[local], beta);
#pragma omp single
      {
        try {
          for (int beta = 0; beta < n_beta_; ++beta) {
            if (!support) {
              read_column(beta, halo, source_column.topRows(halo.size()));
              for (int block = 0; block < block_width; ++block)
                input.col(static_cast<Eigen::Index>(block) * n_beta_ + beta)
                    .head(halo.size()) = source_column.col(block).head(halo.size());
              continue;
            }
            column_rows.clear();
            for (int local = 0; local < static_cast<int>(halo.size()); ++local) {
              if (allowed[static_cast<std::size_t>(beta) * maximum_halo + local])
                column_rows.push_back(halo[local]);
            }
            if (!column_rows.empty())
              read_column(beta, column_rows, source_column.topRows(column_rows.size()));
            for (int block = 0; block < block_width; ++block) {
              auto column = input.col(static_cast<Eigen::Index>(block) * n_beta_ + beta);
              column.head(halo.size()).setZero();
              for (int row = 0; row < static_cast<int>(column_rows.size()); ++row)
                column[halo_index[column_rows[row]]] = source_column(row, block);
            }
          }
          xmvb::profiling::record_matrix_shape(xmvb::profiling::Phase::CiSigma,
              active_rows, team_size);
        } catch (...) { failure = std::current_exception(); }
      }
      if (failure) break;
      beta_reduction.topRows(active_rows).setZero();

#pragma omp for schedule(static)
      for (int beta_target = 0; beta_target < n_beta_; ++beta_target) {
        active_alpha.clear();
        for (int local = 0; local < alpha_count; ++local)
          if (!support || support->contains(alpha_begin + local, beta_target, 1))
            active_alpha.push_back(local);
        const int n_active = static_cast<int>(active_alpha.size());
        if (n_active == 0) continue;
        const int kernel_rows = n_active * block_width;
        intermediate.topRows(kernel_rows).setZero();
        const auto add_beta_source = [&](int beta_source, int pair, double sign) {
          for (int block = 0; block < block_width; ++block) {
            const Eigen::Index column =
                static_cast<Eigen::Index>(block) * n_beta_ + beta_source;
            const int row_base = block * n_active;
            if (!support) {
              intermediate.col(pair).segment(row_base, alpha_count).noalias() +=
                  sign * input.col(column).segment(halo_index[alpha_begin], alpha_count);
            } else {
#pragma omp simd
              for (int local = 0; local < n_active; ++local)
                intermediate(row_base + local, pair) += sign *
                    input(halo_index[alpha_begin] + active_alpha[local], column);
            }
          }
        };
        for (const int occupied : beta_graph.occupied[beta_target])
          add_beta_source(beta_target,
              xmvb::core::PairIndex::packed_pair_index(occupied, occupied),
              1.0);
        const DensityConnections& beta_links = beta_graph.singles[beta_target];
        for (std::size_t link = 0; link < beta_links.size(); ++link)
          add_beta_source(beta_links.sources[link], beta_links.pairs[link],
                          beta_links.signs[link]);

        for (int local_alpha = 0; local_alpha < n_active; ++local_alpha) {
          const int alpha_target = alpha_begin + active_alpha[local_alpha];
          const auto add_alpha_source = [&](int alpha_source, int pair, double sign) {
            if (support && !allowed[static_cast<std::size_t>(beta_target) * maximum_halo +
                                    halo_index[alpha_source]]) return;
            for (int block = 0; block < block_width; ++block) {
              intermediate(block * n_active + local_alpha, pair) +=
                  sign * input(halo_index[alpha_source],
                      static_cast<Eigen::Index>(block) * n_beta_ + beta_target);
            }
          };
          for (const int occupied : alpha_->occupied[alpha_target])
            add_alpha_source(alpha_target,
                xmvb::core::PairIndex::packed_pair_index(occupied, occupied),
                1.0);
          const DensityConnections& alpha_links = alpha_->singles[alpha_target];
          for (std::size_t link = 0; link < alpha_links.size(); ++link)
            add_alpha_source(alpha_links.sources[link], alpha_links.pairs[link],
                             alpha_links.signs[link]);
        }

        transformed.topRows(kernel_rows).noalias() =
            intermediate.topRows(kernel_rows) *
            hamiltonian_pair_kernel_.transpose();

        // Alpha adjoint writes only this beta_target column, owned by the
        // current worker. Beta adjoint uses a thread-private reduction.
        for (int local_alpha = 0; local_alpha < n_active; ++local_alpha) {
          const int alpha_target = alpha_begin + active_alpha[local_alpha];
          const auto add_alpha_image = [&](int alpha_source, int pair, double sign) {
            const int row = halo_index[alpha_source];
            if (support && !allowed[static_cast<std::size_t>(beta_target) * maximum_halo + row])
              return;
            for (int block = 0; block < block_width; ++block)
              output(row, static_cast<Eigen::Index>(block) * n_beta_ + beta_target) +=
                  sign * transformed(block * n_active + local_alpha, pair);
          };
          for (const int occupied : alpha_->occupied[alpha_target])
            add_alpha_image(alpha_target,
                xmvb::core::PairIndex::packed_pair_index(occupied, occupied),
                1.0);
          const DensityConnections& alpha_links = alpha_->singles[alpha_target];
          for (std::size_t link = 0; link < alpha_links.size(); ++link)
            add_alpha_image(alpha_links.sources[link], alpha_links.pairs[link],
                            alpha_links.signs[link]);
        }

        const auto add_beta_image = [&](int beta_source, int pair, double sign) {
          const double* source = transformed.col(pair).data();
          double* target = beta_reduction.col(beta_source).data();
          if (support) {
            const auto* present = allowed.data() +
                static_cast<std::size_t>(beta_source) * maximum_halo + halo_index[alpha_begin];
            for (int block = 0; block < block_width; ++block) {
#pragma omp simd
              for (int local = 0; local < n_active; ++local) {
                const int alpha = active_alpha[local];
                if (present[alpha])
                  target[block * alpha_count + alpha] +=
                      sign * source[block * n_active + local];
              }
            }
            return;
          }
#pragma omp simd
          for (int row = 0; row < active_rows; ++row)
            target[row] += sign * source[row];
        };
        for (const int occupied : beta_graph.occupied[beta_target])
          add_beta_image(beta_target,
              xmvb::core::PairIndex::packed_pair_index(occupied, occupied),
              1.0);
        for (std::size_t link = 0; link < beta_links.size(); ++link)
          add_beta_image(beta_links.sources[link], beta_links.pairs[link],
                         beta_links.signs[link]);
      }

#pragma omp for schedule(static)
      for (int beta = 0; beta < n_beta_; ++beta) {
        double* sum = partial_sigma[0] + static_cast<Eigen::Index>(beta) * maximum_rows;
        for (int source_thread = 1; source_thread < team_size; ++source_thread) {
          const double* source = partial_sigma[source_thread] +
              static_cast<Eigen::Index>(beta) * maximum_rows;
#pragma omp simd
          for (int row = 0; row < active_rows; ++row) sum[row] += source[row];
        }
        for (int block = 0; block < block_width; ++block)
          output.col(static_cast<Eigen::Index>(block) * n_beta_ + beta)
              .segment(halo_index[alpha_begin], alpha_count) +=
              Eigen::Map<const Eigen::VectorXd>(sum + block * alpha_count, alpha_count);
      }

      // The store remains serial and additive. The implicit single barrier
      // prevents the next tile from overwriting a still-consumed output.
#pragma omp single
      {
        try {
          for (int beta = 0; beta < n_beta_; ++beta) {
            if (!support) {
              for (int block = 0; block < block_width; ++block)
                output_column.col(block).head(halo.size()) =
                    output.col(static_cast<Eigen::Index>(block) * n_beta_ + beta).head(halo.size());
              add_column(beta, halo, output_column.topRows(halo.size()));
              continue;
            }
            column_rows.clear();
            for (int local = 0; local < static_cast<int>(halo.size()); ++local)
              if (allowed[static_cast<std::size_t>(beta) * maximum_halo + local])
                column_rows.push_back(halo[local]);
            if (column_rows.empty()) continue;
            for (int block = 0; block < block_width; ++block)
              for (int row = 0; row < static_cast<int>(column_rows.size()); ++row)
                output_column(row, block) = output(halo_index[column_rows[row]],
                    static_cast<Eigen::Index>(block) * n_beta_ + beta);
            add_column(beta, column_rows, output_column.topRows(column_rows.size()));
          }
        } catch (...) { failure = std::current_exception(); }
      }
      if (failure) break;
    }
  }
  if (failure) std::rethrow_exception(failure);
}

}  // namespace xmvb::mcscf
