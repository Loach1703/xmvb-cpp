#include "mcscf/ci/sigma.hpp"

#include <algorithm>
#include <exception>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>

#include "core/openmp.hpp"
#include "core/pair_index.hpp"
#include "core/timing.hpp"
#include "mcscf/ci/product_support.hpp"

namespace xmvb::mcscf {

void DirectCiSigmaAction::apply_tiled(
    int width,
    const std::function<void(const std::vector<int>&, const std::vector<int>&, Eigen::Ref<Eigen::MatrixXd>)>& read_columns,
    const std::function<void(const std::vector<int>&, const std::vector<int>&, const DirectCiColumnAccumulator&)>& add_columns,
    const ProductSupport* support) const {
  profiling::ScopedPhase phase(profiling::Phase::CiSigma);
  if (width <= 0 || !read_columns || !add_columns)
    throw std::invalid_argument("invalid tiled direct-CI callbacks or width");

  const auto& beta_graph = beta_connections();
  const int pairs = hamiltonian_pair_kernel_.rows();
  constexpr std::size_t workspace_bytes = 64ULL * 1024ULL * 1024ULL;
  // Reserve half the slab target for the two channel matrices. The exact
  // halo/CSR payload check below can shorten a batch further. No fixed
  // product count constrains how many beta columns can share a contraction.
  const int row_limit = static_cast<int>(std::max<std::size_t>(1,
      std::min<std::size_t>(std::numeric_limits<int>::max() / width,
          workspace_bytes / (4ULL * sizeof(double) * width * pairs))));
  const int workers = coarse_openmp_thread_count(effective_openmp_thread_count(),
      static_cast<std::size_t>(row_limit) * width * pairs);

  struct Group { int beta, begin, end, capacity; };
  struct Incoming { int group, pair, sign; };
  std::vector<Group> groups;
  std::vector<int> active, group_index, halo, columns, distinct_alpha;
  std::vector<int> halo_index(n_alpha_, -1), column_index(n_beta_, -1);
  std::vector<unsigned char> alpha_seen(n_alpha_, 0);
  std::vector<int> same_group, offsets, positions, target_index;
  std::vector<Incoming> incoming;
  active.reserve(row_limit);
  group_index.reserve(row_limit);
  std::unique_ptr<double[]> workspace;
  std::size_t capacity = 0, panel_size = 0, channel_size = 0;
  int count = 0, n_halo = 0, n_columns = 0;

  // Project the prefix tree onto beta first. Only O(U_beta) column indices
  // survive here; alpha products are generated independently inside a batch.
  std::vector<int> beta_targets;
  if (support) {
    auto cursor = support->cursor(1);
    int beta;
    while (cursor.next_column(beta)) beta_targets.push_back(beta);
  } else {
    for (int beta = 0; beta < n_beta_; ++beta) beta_targets.push_back(beta);
  }
  std::size_t next_beta = 0;
  std::optional<ProductSupport::Cursor> continuation;
  int full_alpha_offset = 0;
  bool continued = false, finished = false;
  std::size_t maximum_neighbors = 1;
  for (const auto& links : alpha_->singles)
    maximum_neighbors = std::max(maximum_neighbors, links.size() + 1);
  // Retain the former two-panel bound after removing the output slab, so
  // batches remain comparable under the same 64 MiB target.
  const auto payload = [&](std::size_t rows, std::size_t h, std::size_t b, std::size_t g) {
    return 2.0L * sizeof(double) * width *
        (static_cast<long double>(h) * b + static_cast<long double>(rows) * pairs) +
        sizeof(int) * static_cast<long double>(h) * g;
  };
  std::exception_ptr failure;
  const bool timing = profiling::phase_totals().enabled;
  using Clock = std::chrono::steady_clock;
  const auto tick = [&]() { return timing ? Clock::now() : Clock::time_point{}; };
  const auto seconds = [](Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double>(end - begin).count();
  };
  double prepare_time = 0, read_time = 0, compute_time = 0, write_time = 0;
  std::size_t batches = 0, source_columns = 0, source_elements = 0, peak_payload = 0;
  Clock::time_point compute_start, prepare_start;
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
          prepare_start = tick();
          for (const int a : halo) halo_index[a] = -1;
          for (const int b : columns) column_index[b] = -1;
          for (const int a : distinct_alpha) alpha_seen[a] = 0;
          distinct_alpha.clear();
          active.clear();
          group_index.clear();
          groups.clear();
          halo.clear();
          columns.clear();
          int reserved = 0;
          while (next_beta + groups.size() < beta_targets.size() && reserved < row_limit) {
            const int beta = beta_targets[next_beta + groups.size()];
            std::size_t b = columns.size();
            b += column_index[beta] < 0;
            for (const int source : beta_graph.singles[beta].sources)
              b += column_index[source] < 0;
            const auto bound = [&](int rows) {
              const auto h = std::min<std::size_t>(n_alpha_,
                  static_cast<std::size_t>(rows) * maximum_neighbors);
              return payload(rows, h, b, groups.size() + 1);
            };
            int slots = std::min(n_alpha_, row_limit - reserved);
            if (!groups.empty() && (slots < n_alpha_ ||
                bound(reserved + slots) > workspace_bytes)) break;
            // Split a single large alpha column without retaining its full
            // list. The same cursor resumes next batch; nothing is regenerated.
            if (groups.empty()) {
              int low = 1, high = slots;
              while (low < high) {
                const int middle = low + (high - low + 1) / 2;
                if (bound(middle) <= workspace_bytes) low = middle;
                else high = middle - 1;
              }
              slots = low;
            }
            const auto include = [](int i, std::vector<int>& ids, std::vector<int>& index) {
              if (index[i] < 0) { index[i] = 0; ids.push_back(i); }
            };
            groups.push_back({beta, reserved, reserved, slots});
            reserved += slots;
            include(beta, columns, column_index);
            for (const int source : beta_graph.singles[beta].sources)
              include(source, columns, column_index);
            if (slots < n_alpha_) break;
          }
          active.resize(reserved);
          continued = false;
          finished = groups.empty();
        } catch (...) { failure = std::current_exception(); }
      }
      if (failure || finished) break;
#pragma omp for schedule(static)
      for (int g = 0; g < static_cast<int>(groups.size()); ++g) {
        try {
          auto& group = groups[g];
          if (support) {
            auto cursor = g == 0 && continuation ? std::move(*continuation) :
                support->column(group.beta, 1);
            int a, b;
            bool label;
            while (group.end - group.begin < group.capacity && cursor.next(a, b, label))
              active[group.end++] = a;
            if (group.capacity < n_alpha_) {
              continued = group.end - group.begin == group.capacity;
              if (continued) continuation.emplace(std::move(cursor));
              else continuation.reset();
            }
          } else {
            const int first_alpha = g == 0 ? full_alpha_offset : 0;
            for (int a = first_alpha; a < n_alpha_ && group.end - group.begin < group.capacity; ++a)
              active[group.end++] = a;
            if (group.capacity < n_alpha_) {
              full_alpha_offset = first_alpha + group.end - group.begin;
              continued = full_alpha_offset < n_alpha_;
              if (!continued) full_alpha_offset = 0;
            }
          }
        } catch (...) {
#pragma omp critical(ci_tile_failure)
          { failure = std::current_exception(); }
        }
      }
#pragma omp single
      {
        try {
          if (!failure) {
            next_beta += groups.size() - static_cast<std::size_t>(continued);
            int packed = 0;
            for (int g = 0; g < static_cast<int>(groups.size()); ++g) {
              auto& group = groups[g];
              const int length = group.end - group.begin;
              for (int i = 0; i < length; ++i)
                active[packed + i] = active[group.begin + i];
              group.begin = packed;
              packed += length;
              group.end = packed;
            }
            active.resize(packed);
            group_index.resize(packed);
            for (int g = 0; g < static_cast<int>(groups.size()); ++g)
              std::fill(group_index.begin() + groups[g].begin,
                        group_index.begin() + groups[g].end, g);
            const auto include = [&](int a) {
              if (halo_index[a] < 0) { halo_index[a] = 0; halo.push_back(a); }
            };
            for (const int a : active) {
              if (alpha_seen[a]) continue;
              alpha_seen[a] = 1;
              distinct_alpha.push_back(a);
              include(a);
              for (const int source : alpha_->singles[a].sources) include(source);
            }
            if (!active.empty()) {
              count = active.size();
              std::sort(halo.begin(), halo.end());
              std::sort(columns.begin(), columns.end());
              n_halo = halo.size();
              n_columns = columns.size();
              for (int i = 0; i < n_halo; ++i) halo_index[halo[i]] = i;
              for (int i = 0; i < n_columns; ++i) column_index[columns[i]] = i;
              panel_size = static_cast<std::size_t>(n_halo) * n_columns * width;
              channel_size = static_cast<std::size_t>(count) * width * pairs;
              const std::size_t required = panel_size + 2 * channel_size;
              if (required > capacity) {
                workspace.reset();
                // Reserve the planned slab once, without value-initializing it.
                // Successively larger batches must not churn large allocations.
                capacity = std::max(required, workspace_bytes / sizeof(double));
                workspace.reset(new double[capacity]);
              }
              Eigen::Map<Eigen::MatrixXd> input(workspace.get(), n_halo,
                  static_cast<Eigen::Index>(n_columns) * width);
              // Each distinct source column is read once per batch. Callbacks
              // mutate the panel LRU and therefore remain outside worker loops.
              const auto read_start = tick();
              read_columns(columns, halo, input);
              const double batch_read_time = seconds(read_start, tick());
              read_time += batch_read_time;
              ++batches;
              source_columns += n_columns;
              source_elements += panel_size;
              peak_payload = std::max(peak_payload, required * sizeof(double) +
                  sizeof(int) * n_halo * groups.size());

              same_group.assign(n_columns, -1);
              target_index.assign(static_cast<std::size_t>(n_halo) * groups.size(), -1);
              offsets.assign(n_columns + 1, 0);
              for (int g = 0; g < static_cast<int>(groups.size()); ++g) {
                same_group[column_index[groups[g].beta]] = g;
                for (int i = groups[g].begin; i < groups[g].end; ++i)
                  target_index[static_cast<std::size_t>(g) * n_halo + halo_index[active[i]]] = i;
                for (const int source : beta_graph.singles[groups[g].beta].sources)
                  ++offsets[column_index[source] + 1];
              }
              for (int i = 1; i <= n_columns; ++i) offsets[i] += offsets[i - 1];
              positions = offsets;
              incoming.resize(offsets.back());
              for (int g = 0; g < static_cast<int>(groups.size()); ++g) {
                const auto& links = beta_graph.singles[groups[g].beta];
                for (int link = 0; link < static_cast<int>(links.size()); ++link)
                  incoming[positions[column_index[links.sources[link]]]++] =
                      {g, links.pairs[link], links.signs[link]};
              }
              profiling::record_matrix_shape(profiling::Phase::CiSigma, count * width, team);
              prepare_time -= batch_read_time;
            }
            compute_start = tick();
            prepare_time += seconds(prepare_start, compute_start);
          }
        } catch (...) { failure = std::current_exception(); }
      }
      if (failure || finished) break;
      if (active.empty()) continue;

      Eigen::Map<const Eigen::MatrixXd> input(workspace.get(), n_halo,
          static_cast<Eigen::Index>(n_columns) * width);
      // Scatter and alpha adjoint traverse channels of a target node. Store
      // that channel vector contiguously, not one cache-line/page per entry.
      using ChannelMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
      Eigen::Map<ChannelMatrix> scattered(workspace.get() + panel_size, count * width, pairs);
      Eigen::Map<ChannelMatrix> transformed(workspace.get() + panel_size + channel_size,
          count * width, pairs);

      const int first = count * width * thread / team;
      const int last = count * width * (thread + 1) / team;
      for (int row = first; row < last; ++row) {
        scattered.row(row).setZero();
        const int local = row % count, v = row / count;
        const int a = active[local], beta = groups[group_index[local]].beta;
        const int column = column_index[beta] * width + v;
        const double value = input(halo_index[a], column);
        for (const int p : alpha_->occupied[a])
          scattered(row, core::PairIndex::packed_pair_index(p, p)) += value;
        for (const int p : beta_graph.occupied[beta])
          scattered(row, core::PairIndex::packed_pair_index(p, p)) += value;
        const auto& alpha_links = alpha_->singles[a];
        for (int link = 0; link < static_cast<int>(alpha_links.size()); ++link)
          scattered(row, alpha_links.pairs[link]) += alpha_links.signs[link] *
              input(halo_index[alpha_links.sources[link]], column);
        const auto& beta_links = beta_graph.singles[beta];
        for (int link = 0; link < static_cast<int>(beta_links.size()); ++link)
          scattered(row, beta_links.pairs[link]) += beta_links.signs[link] *
              input(halo_index[a], column_index[beta_links.sources[link]] * width + v);
      }
      if (last > first)
        transformed.middleRows(first, last - first).noalias() =
            scattered.middleRows(first, last - first) * hamiltonian_pair_kernel_.transpose();
#pragma omp barrier

#pragma omp single
      {
        try {
          // Evaluate E^T K E directly in the caller's output panel. The callback
          // owns panel lifetimes and schedules disjoint beta/RHS/alpha ranges.
          const DirectCiColumnAccumulator accumulate = [&](int source, int v,
              int first, int last, double* destination) {
            const int same = same_group[source];
            for (int local = first; local < last; ++local) {
              const int a = halo[local];
              double sum = 0;
              if (same >= 0) {
                const int* index = target_index.data() + static_cast<std::size_t>(same) * n_halo;
                if (index[local] >= 0) {
                  const int row = v * count + index[local];
                  for (const int p : alpha_->occupied[a])
                    sum += transformed(row, core::PairIndex::packed_pair_index(p, p));
                  for (const int p : beta_graph.occupied[columns[source]])
                    sum += transformed(row, core::PairIndex::packed_pair_index(p, p));
                }
                const auto& links = alpha_->singles[a];
                for (int link = 0; link < static_cast<int>(links.size()); ++link) {
                  const int h = halo_index[links.sources[link]];
                  const int target = h < 0 ? -1 : index[h];
                  if (target >= 0)
                    sum += links.signs[link] * transformed(v * count + target, links.pairs[link]);
                }
              }
              for (int edge = offsets[source]; edge < offsets[source + 1]; ++edge) {
                const auto& link = incoming[edge];
                const int target = target_index[static_cast<std::size_t>(link.group) * n_halo + local];
                if (target >= 0) sum += link.sign * transformed(v * count + target, link.pair);
              }
              destination[a] += sum;
            }
          };
          const auto write_start = tick();
          compute_time += seconds(compute_start, write_start);
          add_columns(columns, halo, accumulate);
          write_time += seconds(write_start, tick());
        } catch (...) { failure = std::current_exception(); }
      }
      if (failure) break;
    }
  }
  if (failure) std::rethrow_exception(failure);
  if (timing)
    std::fprintf(stderr,
        "XMVB_TILE,batches=%zu,source_columns=%zu,source_elements=%zu,peak_payload=%zu,"
        "prepare=%.9f,read=%.9f,compute=%.9f,write=%.9f\n",
        batches, source_columns, source_elements, peak_payload,
        prepare_time, read_time, compute_time, write_time);
}

} // namespace xmvb::mcscf
