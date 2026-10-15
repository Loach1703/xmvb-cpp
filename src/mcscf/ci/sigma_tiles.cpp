#include "mcscf/ci/sigma.hpp"

#include <algorithm>
#include <exception>
#include <limits>
#include <optional>
#include <stdexcept>

#include "core/openmp.hpp"
#include "core/pair_index.hpp"
#include "core/timing.hpp"
#include "mcscf/ci/product_support.hpp"

namespace xmvb::mcscf {

void DirectCiSigmaAction::apply_tiled(
    int width,
    const std::function<void(int, const std::vector<int>&, Eigen::Ref<Eigen::MatrixXd>)>& read_column,
    const std::function<void(int, const std::vector<int>&, const Eigen::Ref<const Eigen::MatrixXd>&)>& add_column,
    const ProductSupport* support) const {
  profiling::ScopedPhase phase(profiling::Phase::CiSigma);
  if (width <= 0 || !read_column || !add_column ||
      static_cast<Eigen::Index>(width) > std::numeric_limits<Eigen::Index>::max() / n_alpha_)
    throw std::invalid_argument("invalid tiled direct-CI callbacks or width");

  const auto& beta_graph = beta_connections();
  const int pairs = hamiltonian_pair_kernel_.rows();
  const int alpha_degree = alpha_->singles.front().size();
  const int beta_degree = beta_graph.singles.front().size();
  constexpr std::size_t workspace_bytes = 64ULL * 1024ULL * 1024ULL;
  // Only the target column's excitation halo and its actual beta neighbors
  // are resident. Nothing here has a dimension proportional to n_beta_.
  int extent = std::min(n_alpha_, 4096);
  const auto halo_bound = [&](int rows) {
    return static_cast<int>(std::min<std::size_t>(n_alpha_,
        static_cast<std::size_t>(rows) * (alpha_degree + 1)));
  };
  const auto payload = [&](int rows) -> long double {
    return sizeof(double) * static_cast<long double>(width) *
        (3.0L * halo_bound(rows) + static_cast<long double>(rows) * (beta_degree + 2 * pairs));
  };
  while (extent > 1 && payload(extent) > workspace_bytes) extent = std::max(1, extent / 2);
  const int max_halo = halo_bound(extent);
  const int workers = coarse_openmp_thread_count(effective_openmp_thread_count(),
      static_cast<std::size_t>(extent) * pairs);

  Eigen::MatrixXd center(max_halo, width), image(max_halo, width), slab(max_halo, width);
  Eigen::MatrixXd neighbors(extent * width, beta_degree);
  Eigen::MatrixXd scattered(extent * width, pairs), transformed(extent * width, pairs);
  std::vector<int> active, halo;
  active.reserve(extent);
  halo.reserve(max_halo);
  std::vector<int> active_index(n_alpha_, -1), halo_index(n_alpha_, -1);
  std::optional<ProductSupport::Cursor> paths;
  if (support) paths.emplace(support->cursor(1));
  int next_a = 0, next_b = 0, beta = 0;
  bool pending = false, finished = false;
  const auto next = [&](int& a, int& b) {
    if (paths) {
      bool in_support;
      return paths->next(a, b, in_support);
    }
    if (next_b == n_beta_) return false;
    a = next_a++;
    b = next_b;
    if (next_a == n_alpha_) { next_a = 0; ++next_b; }
    return true;
  };
  int saved_a = 0, saved_b = 0;
  int count = 0;
  std::exception_ptr failure;
#pragma omp parallel if(workers > 1) num_threads(workers)
  {
    int thread = 0, team = 1;
#ifdef _OPENMP
    thread = omp_get_thread_num();
    team = omp_get_num_threads();
#endif
    for (;;) {
#pragma omp single
      {
        try {
          for (const int a : active) active_index[a] = -1;
          for (const int a : halo) halo_index[a] = -1;
          active.clear();
          halo.clear();
          if (!pending) pending = next(saved_a, saved_b);
          if (!pending) finished = true;
          else {
            beta = saved_b;
            do {
              active.push_back(saved_a);
              pending = next(saved_a, saved_b);
            } while (pending && saved_b == beta && static_cast<int>(active.size()) < extent);
            std::sort(active.begin(), active.end());
            count = active.size();
            const auto include = [&](int a) {
              if (halo_index[a] < 0) { halo_index[a] = 0; halo.push_back(a); }
            };
            for (int i = 0; i < count; ++i) {
              const int a = active[i];
              active_index[a] = i;
              include(a);
              for (const int source : alpha_->singles[a].sources) include(source);
            }
            std::sort(halo.begin(), halo.end());
            for (int i = 0; i < static_cast<int>(halo.size()); ++i) halo_index[halo[i]] = i;
            read_column(beta, halo, center.topRows(halo.size()));
            const auto& links = beta_graph.singles[beta];
            for (int link = 0; link < static_cast<int>(links.size()); ++link) {
              read_column(links.sources[link], active, slab.topRows(count));
              for (int v = 0; v < width; ++v)
                neighbors.col(link).segment(v * count, count) = slab.col(v).head(count);
            }
            profiling::record_matrix_shape(profiling::Phase::CiSigma, count * width, team);
          }
        } catch (...) { failure = std::current_exception(); }
      }
      if (failure || finished) break;

      // Each worker owns complete target rows in both scatter and GEMM.
      const int first = count * width * thread / team;
      const int last = count * width * (thread + 1) / team;
      for (int row = first; row < last; ++row) {
        scattered.row(row).setZero();
        const int local = row % count, v = row / count;
        const int a = active[local];
        const double value = center(halo_index[a], v);
        for (const int p : alpha_->occupied[a])
          scattered(row, core::PairIndex::packed_pair_index(p, p)) += value;
        for (const int p : beta_graph.occupied[beta])
          scattered(row, core::PairIndex::packed_pair_index(p, p)) += value;
        const auto& alpha_links = alpha_->singles[a];
        for (int link = 0; link < static_cast<int>(alpha_links.size()); ++link)
          scattered(row, alpha_links.pairs[link]) += alpha_links.signs[link] *
              center(halo_index[alpha_links.sources[link]], v);
        const auto& beta_links = beta_graph.singles[beta];
        for (int link = 0; link < static_cast<int>(beta_links.size()); ++link)
          scattered(row, beta_links.pairs[link]) +=
              beta_links.signs[link] * neighbors(row, link);
      }
      if (last > first)
        transformed.middleRows(first, last - first).noalias() =
            scattered.middleRows(first, last - first) * hamiltonian_pair_kernel_.transpose();
#pragma omp barrier

      // Reverse alpha links have the same packed pair and fermion sign.
      // Gather by source row: no atomic writes or worker carrier copies.
#pragma omp for schedule(static)
      for (int source = 0; source < static_cast<int>(halo.size()); ++source) {
        const int a = halo[source], local = active_index[a];
        const auto& links = alpha_->singles[a];
        for (int v = 0; v < width; ++v) {
          double sum = 0.0;
          if (local >= 0) {
            const int row = v * count + local;
            for (const int p : alpha_->occupied[a])
              sum += transformed(row, core::PairIndex::packed_pair_index(p, p));
            for (const int p : beta_graph.occupied[beta])
              sum += transformed(row, core::PairIndex::packed_pair_index(p, p));
          }
          for (int link = 0; link < static_cast<int>(links.size()); ++link) {
            const int target = active_index[links.sources[link]];
            if (target >= 0)
              sum += links.signs[link] * transformed(v * count + target, links.pairs[link]);
          }
          image(source, v) = sum;
        }
      }
#pragma omp single
      {
        try {
          add_column(beta, halo, image.topRows(halo.size()));
          const auto& links = beta_graph.singles[beta];
          for (int link = 0; link < static_cast<int>(links.size()); ++link) {
            for (int v = 0; v < width; ++v)
              slab.col(v).head(count) = links.signs[link] *
                  transformed.col(links.pairs[link]).segment(v * count, count);
            add_column(links.sources[link], active, slab.topRows(count));
          }
        } catch (...) { failure = std::current_exception(); }
      }
      if (failure) break;
    }
  }
  if (failure) std::rethrow_exception(failure);
}

} // namespace xmvb::mcscf
