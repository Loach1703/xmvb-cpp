#include "mcscf/ci/sigma.hpp"

#include <algorithm>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <utility>

#include "parallel/openmp.hpp"
#include "core/pair_index.hpp"
#include "core/timing.hpp"
#include "mcscf/ci/product_support.hpp"
#include "mcscf/ci/sigma_plan.hpp"

namespace xmvb::mcscf {
namespace {

using ChannelMatrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;

/** @brief Exclusive scratch lease; cache one idle slab, never one per thread.
 * Concurrent or recursive actions still own disjoint buffers. Earlier returns
 * release their buffer while another lease is active, so the idle cache cannot
 * add a slab to the concurrent working set. Only pointer ownership is moved.
 */
class ChannelWorkspace {
  struct Storage {
    std::unique_ptr<double[]> data;
    std::size_t capacity=0;
    void swap(Storage& other) noexcept {
      data.swap(other.data);
      std::swap(capacity,other.capacity);
    }
  };
  struct Pool {
    std::mutex mutex;
    Storage idle;
    std::size_t active=0;
  };
  static Pool& pool() { static Pool shared; return shared; }
  Storage storage_;
  std::size_t budget_;
  bool reused_=false;
  std::size_t allocations_=0;

public:
  explicit ChannelWorkspace(std::size_t budget):budget_(budget) {
    auto& shared=pool();
    std::lock_guard lock(shared.mutex);
    ++shared.active;
    storage_.swap(shared.idle);
    if (storage_.capacity>budget_/sizeof(double)) {
      storage_.data.reset();
      storage_.capacity=0;
    }
    reused_=storage_.capacity!=0;
  }
  ~ChannelWorkspace() {
    auto& shared=pool();
    std::lock_guard lock(shared.mutex);
    if (--shared.active==0 && storage_.capacity<=budget_/sizeof(double))
      storage_.swap(shared.idle);
  }
  ChannelWorkspace(const ChannelWorkspace&)=delete;
  ChannelWorkspace& operator=(const ChannelWorkspace&)=delete;

  void reserve(std::size_t required) {
    if (required<=storage_.capacity) return;
    storage_.data.reset();
    storage_.capacity=0;
    const auto capacity=required;
    storage_.data=std::make_unique_for_overwrite<double[]>(capacity);
    storage_.capacity=capacity;
    ++allocations_;
  }
  double* get() const noexcept { return storage_.data.get(); }
  double& operator[](std::size_t index) const noexcept { return storage_.data[index]; }
  bool reused() const noexcept { return reused_; }
  std::size_t allocations() const noexcept { return allocations_; }
  std::size_t bytes() const noexcept { return storage_.capacity*sizeof(double); }
};

/** @brief Row-independent product using disjoint, already consumed halo storage. */
void multiply_in_place(Eigen::Ref<ChannelMatrix> values,
    const Eigen::MatrixXd& kernel,Eigen::Ref<ChannelMatrix> scratch) {
  for (Eigen::Index first=0;first<values.rows();first+=scratch.rows()) {
    const auto count=std::min(scratch.rows(),values.rows()-first);
    auto image=scratch.topRows(count);
    image.noalias()=values.middleRows(first,count)*kernel.transpose();
    values.middleRows(first,count)=image;
  }
}

} // namespace

void DirectCiSigmaAction::apply_tiled(
    int width,
    const DirectCiTileReader& read_columns, const DirectCiTileWriter& add_columns,
    const ProductSupport* support, int exchange_parity,
    const DirectCiColumnRelease& release,std::size_t workspace_bytes) const {
  profiling::ScopedPhase phase(profiling::Phase::CiSigma);
  const DirectCiTileRequest request{read_columns,add_columns};
  contract_tiled(width,{&request,1},{},false,support,exchange_parity,release,workspace_bytes);
}

DirectCiIntegralAdjoint DirectCiSigmaAction::apply_tiled_with_adjoint(
    int width, const DirectCiTileReader& response, const DirectCiTileWriter& add_columns,
    const DirectCiTileReader& state, const ProductSupport* support, int exchange_parity,
    const DirectCiColumnRelease& release,std::size_t workspace_bytes) const {
  profiling::ScopedPhase phase(profiling::Phase::CiResponse);
  const DirectCiTileRequest request{response,add_columns};
  auto results=contract_tiled(width,{&request,1},state,true,support,exchange_parity,release,workspace_bytes);
  return std::move(results.front());
}

DirectCiIntegralAdjoint DirectCiSigmaAction::apply_tiled_direction(
    int width, const DirectCiTileReader& response, const DirectCiTileReader& state,
    const DirectCiSigmaAction& direction, const DirectCiTileWriter& add_columns,
    const DirectCiTileWriter& add_forcing, const ProductSupport* support,
    int exchange_parity,const DirectCiColumnRelease& release,std::size_t workspace_bytes) const {
  profiling::ScopedPhase phase(profiling::Phase::CiResponse);
  const DirectCiTileRequest request{response,add_columns,&direction,add_forcing};
  auto results=contract_tiled(width,{&request,1},state,true,support,exchange_parity,release,workspace_bytes);
  return std::move(results.front());
}

std::vector<DirectCiIntegralAdjoint> DirectCiSigmaAction::apply_tiled_responses(
    int width,std::span<const DirectCiTileRequest> requests,const DirectCiTileReader& state,
    const ProductSupport* support,int exchange_parity,
    const DirectCiColumnRelease& release,std::size_t workspace_bytes) const {
  profiling::ScopedPhase phase(profiling::Phase::CiResponse);
  if (!state) throw std::invalid_argument("temporal responses require a common state");
  return contract_tiled(width,requests,state,true,support,exchange_parity,release,workspace_bytes);
}

std::vector<DirectCiIntegralAdjoint> DirectCiSigmaAction::contract_tiled(
    int width,std::span<const DirectCiTileRequest> requests,const DirectCiTileReader& state,
    bool collect_integrals,const ProductSupport* support,int exchange_parity,
    const DirectCiColumnRelease& release,std::size_t workspace_bytes) const {
  if (width<=0 || requests.empty() || (requests.size()>1 && !state))
    throw std::invalid_argument("invalid tiled direct-CI width or response sequence");
  if (exchange_parity < -1 || exchange_parity > 1 ||
      (exchange_parity != 0 && distinct_beta_))
    throw std::invalid_argument("exchange projection requires identical spin spaces");
  bool has_direction=false;
  for (const auto& request:requests) {
    const auto* direction=request.direction;
    if (!request.response || (!request.image && !collect_integrals) ||
        (direction && !state) || (request.forcing && !direction))
      throw std::invalid_argument("invalid tiled direct-CI response callbacks");
    if (direction && (n_orbitals_ != direction->n_orbitals_ ||
                    n_alpha_ != direction->n_alpha_ || n_beta_ != direction->n_beta_ ||
                    alpha_ != direction->alpha_ ||
                    distinct_beta_ != direction->distinct_beta_ ||
                    hamiltonian_pair_kernel_.rows() !=
                        direction->hamiltonian_pair_kernel_.rows()))
      throw std::invalid_argument("direction requires shared direct-CI spin links");
    has_direction=has_direction || direction;
  }

  const auto& beta_graph = beta_connections();
  const int pairs = hamiltonian_pair_kernel_.rows();
  std::vector<Eigen::MatrixXd> kernel_adjoints(requests.size());
  if (collect_integrals) for (auto& gradient:kernel_adjoints)
    gradient=Eigen::MatrixXd::Zero(pairs,pairs);
  // Only these channel slabs coexist: response/state and, for a direction,
  // its output. Response and state refill the same input halo sequentially.
  const bool separate_output=has_direction || requests.size()>1;
  // Earlier CI outputs are consumed before the final orbital request. Only
  // that last request may destroy B to produce its independent forcing image.
  const bool temporal_in_place=requests.size()>1 && std::none_of(
      requests.begin(),requests.end()-1,[](const auto& r) { return r.direction!=nullptr; });
  const std::size_t slab_factor = temporal_in_place ? 2ULL:(separate_output ? 3ULL:2ULL);
  const int row_limit = static_cast<int>(std::max<std::size_t>(1,
      std::min<std::size_t>(std::numeric_limits<int>::max() / width,
          workspace_bytes / (slab_factor * sizeof(double) * width * pairs))));
  const int workers = coarse_openmp_thread_count(effective_openmp_thread_count(),
      static_cast<std::size_t>(row_limit) * width * pairs);

  struct Incoming { int group, pair; double sign; };
  detail::SigmaTilePlan plan(*alpha_,beta_graph,support,exchange_parity,
      {width,pairs,static_cast<int>(slab_factor),workers,temporal_in_place,workspace_bytes},
      static_cast<bool>(release));
  auto& groups=plan.groups;
  auto& active=plan.active;
  auto& group_index=plan.group_index;
  auto& halo=plan.halo;
  auto& columns=plan.columns;
  auto& halo_index=plan.halo_index;
  auto& column_index=plan.column_index;
  auto& active_support=plan.active_support;
  auto& released_columns=plan.released;
  std::vector<int> diagonal_rows,same_group,offsets,positions,target_index;
  std::vector<Eigen::Index> beta_offsets;
  std::vector<Incoming> incoming;
  ChannelWorkspace workspace(workspace_bytes);
  std::size_t panel_size=0,channel_size=0,scratch_size=0;
  int count=0,n_halo=0,n_columns=0;
  bool finished=false;
  std::exception_ptr failure;
  const bool timing = profiling::phase_totals().enabled;
  using Clock = std::chrono::steady_clock;
  const auto tick = [&]() { return timing ? Clock::now() : Clock::time_point{}; };
  const auto seconds = [](Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double>(end - begin).count();
  };
  double prepare_time = 0, read_time = 0, compute_time = 0, write_time = 0, gram_time = 0;
  std::size_t gram_tiles=0;
  std::size_t gram_scratch_bytes=0;
  int gram_partitions=0,gram_extent=0;
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
          prepare_start=tick();
          finished=!plan.next();
          diagonal_rows.clear();
          if (!failure) {
            for (int i=0;i<static_cast<int>(active.size());++i) {
              unsupported_targets+=!active_support[i];
              if (exchange_parity && active[i]==groups[group_index[i]].beta)
                diagonal_rows.push_back(i);
            }
            if (!active.empty()) {
              count = active.size();
              n_halo = halo.size();
              n_columns = columns.size();
              panel_size = static_cast<std::size_t>(n_halo) * n_columns * width;
              channel_size = static_cast<std::size_t>(count) * width * pairs;
              scratch_size=temporal_in_place ? std::max(panel_size,
                  static_cast<std::size_t>(std::min(team,count*width))*pairs):panel_size;
              const std::size_t required = scratch_size + slab_factor * channel_size;
              workspace.reserve(required);
              Eigen::Map<Eigen::MatrixXd> input(workspace.get(), n_halo,
                  static_cast<Eigen::Index>(n_columns) * width);
              // Each distinct source column is read once per batch. Callbacks
              // mutate the panel LRU and therefore remain outside worker loops.
              const auto read_start = tick();
              if (timing && batches==0)
                std::fprintf(stderr,"XMVB_TILE_PROGRESS,event=first_read,rows=%d,columns=%d,halo=%d\n",
                    count,n_columns,n_halo);
              requests.front().response(columns, halo, input);
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
              beta_map_bytes=std::max(beta_map_bytes,beta_offsets.capacity()*sizeof(Eigen::Index));
              for (int g = 0; g < static_cast<int>(groups.size()); ++g) {
                groups[g].column=column_index[groups[g].beta]*width;
                const auto& links = beta_graph.singles(groups[g].beta);
                for (int link = 0; link < static_cast<int>(links.size()); ++link) {
                  incoming[positions[column_index[links[link].source]]++] =
                      {g, links[link].pair, links[link].sign};
                  // The same beta edge serves every alpha target and RHS.
                  // Offsets follow the shared runtime budget, not a 32-bit
                  // assumption tied to the former 64 MiB halo.
                  beta_offsets[g*beta_graph.n_singles()+link]=static_cast<Eigen::Index>(
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

      Eigen::Map<const Eigen::MatrixXd> input(workspace.get(), n_halo,
          static_cast<Eigen::Index>(n_columns) * width);
      // Keep each node's channel vector contiguous through scatter and GEMM.
      Eigen::Map<ChannelMatrix> scattered(workspace.get() + scratch_size, count * width, pairs);
      Eigen::Map<ChannelMatrix> state_channels(workspace.get() + scratch_size + channel_size,
          count * width, pairs);
      Eigen::Map<ChannelMatrix> transformed(workspace.get() + scratch_size +
          (temporal_in_place ? 0:(separate_output ? 2:1)) * channel_size,
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
      // Orbit weights are 2 except on a=b. Apply them to the small Gram,
      // not to B: temporal responses and deltaK*B share its unweighted slab.
      const auto subtract_diagonal=[&](auto&& gram,int begin,int end,int row,int col) {
        for (int v=0;v<width;++v) {
          const int lo=std::max(0,begin-v*count),hi=std::min(count,end-v*count);
          if (lo>=hi) continue;
          auto entry=std::lower_bound(diagonal_rows.begin(),diagonal_rows.end(),lo);
          for (;entry!=diagonal_rows.end() && *entry<hi;++entry) {
            const int r=v*count+*entry;
            gram.noalias()-=scattered.row(r).segment(row,gram.rows()).transpose()*
                state_channels.row(r).segment(col,gram.cols());
          }
        }
      };
      for (std::size_t request_id=0;request_id<requests.size();++request_id) {
        const auto& request=requests[request_id];
        const auto& add_columns=request.image;
        const auto& add_forcing=request.forcing;
        const auto* direction=request.direction;
        auto& kernel_adjoint=kernel_adjoints[request_id];
        if (request_id>0) {
#pragma omp single
          {
            try {
              const auto begin=tick();
              Eigen::Map<Eigen::MatrixXd> values(workspace.get(),n_halo,n_columns*width);
              request.response(columns,halo,values);
              read_time+=seconds(begin,tick());
              source_columns+=n_columns;
              source_elements+=panel_size;
              compute_start=tick();
            } catch (...) { failure=std::current_exception(); }
          }
          if (failure) break;
        }
        for (int row = first; row < last; ++row) scatter(row,scattered);
        if (collect_integrals) {
          if (state && request_id==0) {
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
          } else if (!state && last > first) {
            state_channels.middleRows(first,last-first) = scattered.middleRows(first,last-first);
          }
          // B remains unweighted and immutable for every temporal response.
#pragma omp barrier
          // G=sum_t A_t^T B_t. The consumed input halo can hold the local Gram
          // images without allocating or copying channel slabs. If it cannot
          // hold full local images, own output blocks directly instead: tiny
          // halos must not launch one barrier per scalar Gram entry.
          const int available=std::min(team,count*width);
          const auto full_area=static_cast<std::size_t>(pairs)*pairs;
          const bool row_partition=available>1 && full_area<=scratch_size/available;
          const int parts=row_partition ? available:1;
          int axis=1;
          if (!row_partition) while (axis*axis<team) axis*=2;
          const int extent=(pairs+axis-1)/axis;
          const int blocks=(pairs+extent-1)/extent;
          const auto tile_area=static_cast<std::size_t>(extent)*extent;
          if (timing) {
#pragma omp single
            {
              gram_start=tick();
              gram_partitions=std::max(gram_partitions,parts);
              gram_extent=std::max(gram_extent,extent);
              if (row_partition) gram_scratch_bytes=std::max(gram_scratch_bytes,
                  parts*tile_area*sizeof(double));
            }
          }
          if (row_partition) {
            if (thread<parts) {
              const int begin=count*width*thread/parts,end=count*width*(thread+1)/parts;
              Eigen::Map<Eigen::MatrixXd> partial(
                  workspace.get()+thread*tile_area,pairs,pairs);
              partial.noalias()=(exchange_parity ? 2.0:1.0)*
                  scattered.middleRows(begin,end-begin).transpose()*
                  state_channels.middleRows(begin,end-begin);
              if (exchange_parity) subtract_diagonal(partial,begin,end,0,0);
            }
#pragma omp barrier
#pragma omp for schedule(static)
            for (int entry=0;entry<pairs*pairs;++entry) {
              double value=0;
              for (int part=0;part<parts;++part)
                value+=workspace[part*tile_area+entry];
              kernel_adjoint.data()[entry]+=value;
            }
          } else {
#pragma omp for schedule(static)
            for (int block=0;block<blocks*blocks;++block) {
              const int row=(block/blocks)*extent,col=(block%blocks)*extent;
              const int rows=std::min(extent,pairs-row),cols=std::min(extent,pairs-col);
              auto gram=kernel_adjoint.block(row,col,rows,cols);
              gram.noalias()+=(exchange_parity ? 2.0:1.0)*
                  scattered.middleCols(row,rows).transpose()*state_channels.middleCols(col,cols);
              if (exchange_parity) subtract_diagonal(gram,0,count*width,row,col);
            }
          }
          if (timing) {
#pragma omp single
            { gram_time+=seconds(gram_start,tick()); gram_tiles+=parts*blocks*blocks; }
          }
        }
        if (temporal_in_place && last>first) {
          const int owners=std::min(team,count*width);
          const int owner=count*width<team ? first:thread;
          const Eigen::Index rows=scratch_size/(static_cast<std::size_t>(owners)*pairs);
          Eigen::Map<ChannelMatrix> scratch(workspace.get()+
              static_cast<std::size_t>(owner)*rows*pairs,rows,pairs);
          if (add_columns) multiply_in_place(scattered.middleRows(first,last-first),
              hamiltonian_pair_kernel_,scratch);
          if (direction && (add_columns || add_forcing)) {
            multiply_in_place(state_channels.middleRows(first,last-first),
                direction->hamiltonian_pair_kernel_,scratch);
            if (add_columns) transformed.middleRows(first,last-first)+=
                state_channels.middleRows(first,last-first);
          }
        }
        if (!temporal_in_place && add_columns && last > first)
          transformed.middleRows(first, last - first).noalias() =
              scattered.middleRows(first, last - first) * hamiltonian_pair_kernel_.transpose();
        if (!temporal_in_place && direction && (add_columns || add_forcing) && last > first) {
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
            if (active[local] != groups[group_index[local]].beta) {
              if (temporal_in_place) state_channels.row(row)*=2.0;
              else scattered.row(row) *= 2.0;
            }
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
            if (add_forcing) write_image(add_forcing,
                temporal_in_place ? state_channels.data():scattered.data());
            // Retire inside the final writer's existing single/barrier.
            // Testing the shared retirement vector outside this region races
            // the next plan.next(): one worker could enter a conditional
            // barrier after another has already advanced to the next tile.
            if (request_id+1==requests.size() && release && !released_columns.empty())
              release(released_columns);
            write_time += seconds(write_start, tick());
            if (timing && (batches==1 || batches%256==0))
              std::fprintf(stderr,"XMVB_TILE_PROGRESS,batches=%zu,targets=%zu,prepare=%.9f,read=%.9f,compute=%.9f,write=%.9f\n",
                  batches,targets,prepare_time,read_time,compute_time,write_time);
          } catch (...) { failure = std::current_exception(); }
        }
        if (failure) break;
      } // temporal responses share the generated graph and accepted channels
      if (failure) break;
    }
  }
  if (failure) std::rethrow_exception(failure);
  if (timing)
    std::fprintf(stderr,
        "XMVB_TILE,batches=%zu,source_columns=%zu,source_elements=%zu,peak_payload=%zu,"
        "prepare=%.9f,read=%.9f,compute=%.9f,write=%.9f,targets=%zu,exchange=%d,"
        "skipped_diagonal_scatters=%zu,gram=%.9f,gram_tiles=%zu,"
        "identity_halo_targets=%zu,beta_map_bytes=%zu,gram_partitions=%d,"
        "gram_extent=%d,gram_scratch_bytes=%zu,workspace_reused=%d,"
        "workspace_allocations=%zu,workspace_bytes=%zu\n",
        batches, source_columns, source_elements, peak_payload,
        prepare_time, read_time, compute_time, write_time, targets, exchange_parity,
        unsupported_targets * static_cast<std::size_t>(width) *
            (requests.size() + static_cast<std::size_t>(collect_integrals && static_cast<bool>(state))),
        gram_time,gram_tiles,identity_halo_targets,beta_map_bytes,
        gram_partitions,gram_extent,gram_scratch_bytes,
        static_cast<int>(workspace.reused()),workspace.allocations(),workspace.bytes());
  std::vector<DirectCiIntegralAdjoint> results(requests.size());
  if (collect_integrals) for (std::size_t r=0;r<requests.size();++r)
    results[r]=pullback_channel_kernel(kernel_adjoints[r]);
  return results;
}

} // namespace xmvb::mcscf
