#include "mcscf/ci/sigma.hpp"

#include <algorithm>
#include <exception>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>

#include "parallel/openmp.hpp"
#include "core/pair_index.hpp"
#include "core/timing.hpp"
#include "mcscf/ci/product_support.hpp"

namespace xmvb::mcscf {

void DirectCiSigmaAction::apply_tiled(
    int width,
    const DirectCiTileReader& read_columns, const DirectCiTileWriter& add_columns,
    const ProductSupport* support, int exchange_parity) const {
  profiling::ScopedPhase phase(profiling::Phase::CiSigma);
  contract_tiled(width,read_columns,add_columns,{},false,support,exchange_parity);
}

DirectCiIntegralAdjoint DirectCiSigmaAction::apply_tiled_with_adjoint(
    int width, const DirectCiTileReader& response, const DirectCiTileWriter& add_columns,
    const DirectCiTileReader& state, const ProductSupport* support, int exchange_parity) const {
  profiling::ScopedPhase phase(profiling::Phase::CiResponse);
  return contract_tiled(width,response,add_columns,state,true,support,exchange_parity);
}

DirectCiIntegralAdjoint DirectCiSigmaAction::apply_tiled_direction(
    int width, const DirectCiTileReader& response, const DirectCiTileReader& state,
    const DirectCiSigmaAction& direction, const DirectCiTileWriter& add_columns,
    const DirectCiTileWriter& add_forcing, const ProductSupport* support,
    int exchange_parity) const {
  profiling::ScopedPhase phase(profiling::Phase::CiResponse);
  return contract_tiled(width,response,add_columns,state,true,support,exchange_parity,
                        &direction,add_forcing);
}

DirectCiIntegralAdjoint DirectCiSigmaAction::contract_tiled(
    int width, const DirectCiTileReader& read_columns, const DirectCiTileWriter& add_columns,
    const DirectCiTileReader& state, bool collect_integrals,
    const ProductSupport* support, int exchange_parity,
    const DirectCiSigmaAction* direction, const DirectCiTileWriter& add_forcing) const {
  if (width <= 0 || !read_columns || (!add_columns && !collect_integrals) ||
      (direction && !state) || (add_forcing && !direction))
    throw std::invalid_argument("invalid tiled direct-CI callbacks or width");
  if (exchange_parity < -1 || exchange_parity > 1 ||
      (exchange_parity != 0 && distinct_beta_))
    throw std::invalid_argument("exchange projection requires identical spin spaces");
  if (direction && (n_orbitals_ != direction->n_orbitals_ ||
                    n_alpha_ != direction->n_alpha_ || n_beta_ != direction->n_beta_ ||
                    alpha_ != direction->alpha_ ||
                    distinct_beta_ != direction->distinct_beta_ ||
                    hamiltonian_pair_kernel_.rows() !=
                        direction->hamiltonian_pair_kernel_.rows()))
    throw std::invalid_argument("direction requires shared direct-CI spin links");

  const auto& beta_graph = beta_connections();
  const int pairs = hamiltonian_pair_kernel_.rows();
  Eigen::MatrixXd kernel_adjoint = Eigen::MatrixXd::Zero(
      collect_integrals ? pairs : 0,collect_integrals ? pairs : 0);
  constexpr std::size_t workspace_bytes = 64ULL * 1024ULL * 1024ULL;
  // The ordinary path retains its two-slab limit. The direction path reserves
  // room for a third shared channel slab under the same 64 MiB target.
  const std::size_t slab_factor = direction ? 6ULL : 4ULL;
  const int row_limit = static_cast<int>(std::max<std::size_t>(1,
      std::min<std::size_t>(std::numeric_limits<int>::max() / width,
          workspace_bytes / (slab_factor * sizeof(double) * width * pairs))));
  const int workers = coarse_openmp_thread_count(effective_openmp_thread_count(),
      static_cast<std::size_t>(row_limit) * width * pairs);

  struct Group { int beta, begin, end, capacity, limit; int column=0; };
  struct Incoming { int group, pair; double sign; };
  std::vector<Group> groups;
  std::vector<int> active, group_index, halo, columns, distinct_alpha;
  std::vector<int> halo_index(n_alpha_, -1), column_index(n_beta_, -1);
  std::vector<unsigned char> alpha_seen(n_alpha_, 0), active_support;
  std::vector<int> same_group, offsets, positions, target_index,beta_offsets;
  std::vector<Incoming> incoming;
  active.reserve(row_limit);
  active_support.reserve(row_limit);
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
  if (exchange_parity < 0)
    beta_targets.erase(std::remove(beta_targets.begin(), beta_targets.end(), 0),
                       beta_targets.end());
  std::size_t next_beta = 0;
  std::optional<ProductSupport::Cursor> continuation;
  int full_alpha_offset = 0;
  bool continued = false, finished = false;
  const std::size_t maximum_neighbors=alpha_->n_singles()+1;
  // Retain the former two-panel bound for ordinary requests. Direction
  // requests add only the third channel slab to this conservative estimate.
  const auto payload = [&](std::size_t rows, std::size_t h, std::size_t b, std::size_t g) {
    return sizeof(double) * width *
        (2.0L * static_cast<long double>(h) * b +
         static_cast<long double>(direction ? 3 : 2) * rows * pairs) +
        sizeof(int) * static_cast<long double>(h) * g;
  };
  std::exception_ptr failure;
  const bool timing = profiling::phase_totals().enabled;
  using Clock = std::chrono::steady_clock;
  const auto tick = [&]() { return timing ? Clock::now() : Clock::time_point{}; };
  const auto seconds = [](Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double>(end - begin).count();
  };
  double prepare_time = 0, read_time = 0, compute_time = 0, write_time = 0, gram_time = 0;
  std::size_t gram_tiles=0;
  std::size_t batches = 0, source_columns = 0, source_elements = 0, peak_payload = 0;
  std::size_t targets = 0, unsupported_targets = 0, identity_halo_targets=0,beta_map_bytes=0;
  Clock::time_point compute_start, prepare_start, gram_start;
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
            // Reserve only possible exchange representatives, not the absent
            // lower triangle. This reduces batches without increasing memory.
            const int limit = exchange_parity ? beta + (exchange_parity > 0) : n_alpha_;
            std::size_t b = columns.size();
            b += column_index[beta] < 0;
            for (const auto& edge:beta_graph.singles(beta))
              b += column_index[edge.source] < 0;
            const auto bound = [&](int rows) {
              const auto h = std::min<std::size_t>(n_alpha_,
                  static_cast<std::size_t>(rows) * maximum_neighbors);
              return payload(rows, h, b, groups.size() + 1);
            };
            int slots = std::min(limit, row_limit - reserved);
            if (!groups.empty() && (slots < limit ||
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
            groups.push_back({beta, reserved, reserved, slots, limit});
            reserved += slots;
            include(beta, columns, column_index);
            for (const auto& edge:beta_graph.singles(beta))
              include(edge.source, columns, column_index);
            if (slots < limit) break;
          }
          active.resize(reserved);
          active_support.resize(reserved);
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
            while (group.end - group.begin < group.capacity && cursor.next(a, b, label)) {
              if (exchange_parity && (a > b || (a == b && exchange_parity < 0))) continue;
              active[group.end] = a;
              active_support[group.end++] = label;
            }
            if (group.capacity < group.limit) {
              continued = group.end - group.begin == group.capacity;
              if (continued) continuation.emplace(std::move(cursor));
              else continuation.reset();
            }
          } else {
            const int first_alpha = g == 0 ? full_alpha_offset : 0;
            const int alpha_end = exchange_parity ? group.beta + (exchange_parity > 0) : n_alpha_;
            for (int a = first_alpha; a < alpha_end && group.end - group.begin < group.capacity; ++a) {
              active[group.end] = a;
              active_support[group.end++] = 1;
            }
            if (group.capacity < group.limit) {
              full_alpha_offset = first_alpha + group.end - group.begin;
              continued = full_alpha_offset < alpha_end;
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
              for (int i = 0; i < length; ++i) {
                active[packed + i] = active[group.begin + i];
                active_support[packed + i] = active_support[group.begin + i];
                unsupported_targets += !active_support[packed + i];
              }
              group.begin = packed;
              packed += length;
              group.end = packed;
            }
            active.resize(packed);
            active_support.resize(packed);
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
              for (const auto& edge:alpha_->singles(a)) include(edge.source);
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
              const std::size_t required = panel_size +
                  (direction ? 3 : 2) * channel_size;
              if (direction && required * sizeof(double) +
                      sizeof(int) * static_cast<std::size_t>(n_halo) * groups.size() >
                      workspace_bytes)
                throw std::length_error("tiled direction exceeds 64 MiB workspace");
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
              targets += count;
              if (n_halo==n_alpha_) identity_halo_targets+=count;
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
                for (const auto& edge:beta_graph.singles(groups[g].beta))
                  ++offsets[column_index[edge.source] + 1];
              }
              for (int i = 1; i <= n_columns; ++i) offsets[i] += offsets[i - 1];
              positions = offsets;
              incoming.resize(offsets.back());
              beta_offsets.resize(groups.size()*beta_graph.n_singles());
              beta_map_bytes=std::max(beta_map_bytes,beta_offsets.capacity()*sizeof(int));
              for (int g = 0; g < static_cast<int>(groups.size()); ++g) {
                groups[g].column=column_index[groups[g].beta]*width;
                const auto& links = beta_graph.singles(groups[g].beta);
                for (int link = 0; link < static_cast<int>(links.size()); ++link) {
                  incoming[positions[column_index[links[link].source]]++] =
                      {g, links[link].pair, links[link].sign};
                  // The same beta edge serves every alpha target and RHS.
                  // A halo is bounded by 64 MiB, so its double offset fits int.
                  beta_offsets[g*beta_graph.n_singles()+link]=static_cast<int>(
                      static_cast<std::size_t>(n_halo)*column_index[links[link].source]*width);
                }
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
      // Keep each node's channel vector contiguous through scatter and GEMM.
      using ChannelMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
      Eigen::Map<ChannelMatrix> scattered(workspace.get() + panel_size, count * width, pairs);
      Eigen::Map<ChannelMatrix> state_channels(workspace.get() + panel_size + channel_size,
          count * width, pairs);
      Eigen::Map<ChannelMatrix> transformed(workspace.get() + panel_size +
          (direction ? 2 : 1) * channel_size,
          count * width, pairs);

      const int first = count * width * thread / team;
      const int last = count * width * (thread + 1) / team;
      const bool identity_halo=n_halo==n_alpha_;
      const auto scatter = [&](int row, Eigen::Ref<ChannelMatrix> target) {
        target.row(row).setZero();
        const int local = row % count, v = row / count;
        const int a=active[local],g=group_index[local],beta=groups[g].beta;
        const int h=identity_halo ? a:halo_index[a];
        const double* column=input.col(groups[g].column+v).data();
        // A depth-one target outside S has no own coefficient, although its
        // one-body neighbors can still contribute to E*c and must be kept.
        if (active_support[local]) {
          const double value = column[h];
          for (const int p : alpha_->occupied(a))
            target(row, core::PairIndex::packed_pair_index(p, p)) += value;
          for (const int p : beta_graph.occupied(beta))
            target(row, core::PairIndex::packed_pair_index(p, p)) += value;
        }
        const auto& alpha_links = alpha_->singles(a);
        for (int link = 0; link < static_cast<int>(alpha_links.size()); ++link)
          target(row, alpha_links[link].pair) += alpha_links[link].sign *
              column[identity_halo ? alpha_links[link].source:halo_index[alpha_links[link].source]];
        const auto& beta_links = beta_graph.singles(beta);
        const double* source=input.data()+h+static_cast<std::size_t>(n_halo)*v;
        for (int link = 0; link < static_cast<int>(beta_links.size()); ++link)
          target(row, beta_links[link].pair) += beta_links[link].sign *
              source[beta_offsets[g*beta_graph.n_singles()+link]];
      };
      for (int row = first; row < last; ++row) scatter(row,scattered);
      if (collect_integrals) {
        if (state) {
          // The response scatter is complete. Refill the same input halo for
          // the unweighted state, then scatter it into the second channel slab.
#pragma omp barrier
#pragma omp single
          {
            try {
              const auto begin = tick();
              Eigen::Map<Eigen::MatrixXd> values(workspace.get(),n_halo,n_columns*width);
              state(columns,halo,values);
              const double elapsed=seconds(begin,tick());
              read_time += elapsed;
              compute_time -= elapsed;
              source_columns += n_columns;
              source_elements += panel_size;
            } catch (...) { failure = std::current_exception(); }
          }
          if (failure) break;
          for (int row = first; row < last; ++row) scatter(row,state_channels);
        } else if (last > first) {
          state_channels.middleRows(first,last-first) = scattered.middleRows(first,last-first);
        }
        if (exchange_parity)
          for (int row = first; row < last; ++row)
            if (active[row % count] != groups[group_index[row % count]].beta)
              state_channels.row(row) *= 2.0;
#pragma omp barrier
        // G(I,J)=A(:,I)^T B(:,J): a square output grid avoids each narrow
        // row-stripe worker rereading the entire right-channel slab. Views
        // share existing storage; every output block has exactly one owner.
        int axis=1;
        while (axis*axis<team) axis*=2;
        const int extent=(pairs+axis-1)/axis;
        const int blocks=(pairs+extent-1)/extent;
        if (timing) {
#pragma omp single
          { gram_start=tick(); }
        }
#pragma omp for schedule(static)
        for (int block=0;block<blocks*blocks;++block) {
          const int row=(block/blocks)*extent,col=(block%blocks)*extent;
          const int rows=std::min(extent,pairs-row),cols=std::min(extent,pairs-col);
          kernel_adjoint.block(row,col,rows,cols).noalias()+=
              scattered.middleCols(row,rows).transpose()*state_channels.middleCols(col,cols);
        }
        if (timing) {
#pragma omp single
          { gram_time+=seconds(gram_start,tick()); gram_tiles+=blocks*blocks; }
        }
      }
      if (add_columns && last > first)
        transformed.middleRows(first, last - first).noalias() =
            scattered.middleRows(first, last - first) * hamiltonian_pair_kernel_.transpose();
      if (direction && last > first) {
        // The adjoint used orbit multiplicities. The directional image must
        // use the original state; reconstruction is applied only once below.
        if (exchange_parity)
          for (int row = first; row < last; ++row)
            if (active[row % count] != groups[group_index[row % count]].beta)
              state_channels.row(row) *= 0.5;
        // Base K*A has consumed A. Reuse its slab for deltaK*B, then add to C.
        scattered.middleRows(first, last - first).noalias() =
            state_channels.middleRows(first, last - first) *
            direction->hamiltonian_pair_kernel_.transpose();
        if (add_columns)
          transformed.middleRows(first, last - first) +=
              scattered.middleRows(first, last - first);
      }
      // E commutes with spin exchange. In a fixed-parity bilinear form,
      // off-diagonal exchange orbits have two identical contributions.
      if (add_columns && exchange_parity)
        for (int row = first; row < last; ++row) {
          const int local = row % count;
          if (active[local] != groups[group_index[local]].beta)
            transformed.row(row) *= 2.0;
        }
      if (add_forcing && exchange_parity)
        for (int row = first; row < last; ++row) {
          const int local = row % count;
          if (active[local] != groups[group_index[local]].beta)
            scattered.row(row) *= 2.0;
        }
#pragma omp barrier

#pragma omp single
      {
        try {
          // Evaluate E^T K E directly in the caller's output panel. The callback
          // owns panel lifetimes and schedules disjoint beta/RHS/alpha ranges.
          const auto write_image = [&](const DirectCiTileWriter& writer,
                                       const double* channel_values) {
            const DirectCiColumnAccumulator accumulate = [&,channel_values](int source, int v,
              int first,int last,double* destination,const int* rows,int n_rows) {
            if (first==last) return;
            int begin=first,end=last;
            if (rows) {
              begin=std::lower_bound(rows,rows+n_rows,halo[first])-rows;
              end=std::upper_bound(rows+begin,rows+n_rows,halo[last-1])-rows;
              if (begin==end) return;
            }
            const int same = same_group[source];
            if (same >= 0) {
              const int* index = target_index.data() + static_cast<std::size_t>(same) * n_halo;
              for (int i=begin;i<end;++i) {
                const int a=rows ? rows[i]:halo[i];
                const int local=rows ? (identity_halo ? a:halo_index[a]):i;
                if (local<first || local>=last) continue;
                double sum = 0;
                if (index[local] >= 0) {
                  const int row = v * count + index[local];
                  for (const int p : alpha_->occupied(a))
                    sum += channel_values[static_cast<std::size_t>(row) * pairs +
                        core::PairIndex::packed_pair_index(p, p)];
                  for (const int p : beta_graph.occupied(columns[source]))
                    sum += channel_values[static_cast<std::size_t>(row) * pairs +
                        core::PairIndex::packed_pair_index(p, p)];
                }
                const auto& links = alpha_->singles(a);
                for (int link = 0; link < static_cast<int>(links.size()); ++link) {
                  const int h = identity_halo ? links[link].source:halo_index[links[link].source];
                  const int target = h < 0 ? -1 : index[h];
                  if (target >= 0)
                    sum += links[link].sign *
                        channel_values[static_cast<std::size_t>(v * count + target) * pairs +
                            links[link].pair];
                }
                destination[rows ? i:a] += sum;
              }
            }
            // One beta edge fixes the channel, sign and target-row map for
            // the entire alpha range. Do not rediscover them for every pixel.
            for (int edge = offsets[source]; edge < offsets[source + 1]; ++edge) {
              const auto& link = incoming[edge];
              const int* index = target_index.data() + static_cast<std::size_t>(link.group) * n_halo;
              const double* values = channel_values +
                  static_cast<std::size_t>(v) * count * pairs + link.pair;
              const double sign = link.sign;
              for (int i=begin;i<end;++i) {
                const int a=rows ? rows[i]:halo[i];
                const int local=rows ? (identity_halo ? a:halo_index[a]):i;
                if (local<first || local>=last) continue;
                const int target=index[local];
                if (target>=0)
                  destination[rows ? i:a]+=sign*values[static_cast<std::size_t>(target)*pairs];
              }
            }
          };
            writer(columns, halo, accumulate);
          };
          const auto write_start = tick();
          compute_time += seconds(compute_start, write_start);
          if (add_columns) write_image(add_columns, transformed.data());
          if (add_forcing) write_image(add_forcing, scattered.data());
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
        "prepare=%.9f,read=%.9f,compute=%.9f,write=%.9f,targets=%zu,exchange=%d,"
        "skipped_diagonal_scatters=%zu,gram=%.9f,gram_tiles=%zu,"
        "identity_halo_targets=%zu,beta_map_bytes=%zu\n",
        batches, source_columns, source_elements, peak_payload,
        prepare_time, read_time, compute_time, write_time, targets, exchange_parity,
        unsupported_targets * static_cast<std::size_t>(width) *
            (1 + static_cast<std::size_t>(collect_integrals && static_cast<bool>(state))),
        gram_time,gram_tiles,identity_halo_targets,beta_map_bytes);
  return collect_integrals ? pullback_channel_kernel(kernel_adjoint) : DirectCiIntegralAdjoint{};
}

} // namespace xmvb::mcscf
