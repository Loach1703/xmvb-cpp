#pragma once

#include <algorithm>
#include <functional>
#include <initializer_list>
#include <list>
#include <memory>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

#include "parallel/openmp.hpp"
#include "core/timing.hpp"
#include "core/workspace.hpp"

namespace xmvb::vb::detail {

/** @brief Existing structure-action workspace; panels share equal quarters. */
inline constexpr std::size_t kStructureTileWorkspaceBytes=256ULL*1024*1024;

/** @brief Compare exact beta contractions, independently of residency.
 * Row sampling costs touched_products * sampled_beta_columns multiply-adds. Applying the
 * ordered exterior factors costs n_alpha * (n_beta + shear_pairs). Alpha
 * transforms and structure incidence are common to both orders. Neither
 * order assumes closure of the selected structure space.
 */
inline bool prefer_ordered_beta(std::size_t touched_products, int n_alpha,
    int n_beta, std::size_t shear_pairs,int sampled_beta_columns) noexcept {
  return static_cast<long double>(touched_products)*sampled_beta_columns >
      static_cast<long double>(n_alpha)*(n_beta+static_cast<long double>(shear_pairs));
}

/** @brief Exact generated support of beta columns, shared by primal/image tiles. */
struct ColumnPattern {
  bool full=false;
  Eigen::VectorXi offsets;
  Eigen::VectorXi rows;
  int columns() const noexcept { return offsets.size()-1; }
  std::size_t bytes() const noexcept { return sizeof(int)*(offsets.size()+rows.size()); }
};

/** @brief Values only on the symbolic support; RHSs share the row indices. */
struct ColumnTile {
  std::shared_ptr<const ColumnPattern> pattern;
  Eigen::MatrixXd values;

  /** @brief Conservative retained payload, charging even shared metadata. */
  std::size_t bytes() const noexcept { return sizeof(double)*values.size()+pattern->bytes(); }

  /** @brief One temporary full-row panel for the ordered exterior adjoint. */
  Eigen::MatrixXd unpack(int n_alpha) const {
    const int cols=pattern->columns(),width=values.cols();
    if (pattern->full)
      return Eigen::Map<const Eigen::MatrixXd>(values.data(),n_alpha,cols*width);
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(n_alpha,cols*width);
    const int workers=coarse_openmp_thread_count(cooperative_openmp_thread_count(),result.size());
    cooperative_openmp_for(cols*width,workers,[&](std::ptrdiff_t packed) {
      const int b=packed%cols,v=packed/cols;
      for (int i=pattern->offsets[b];i<pattern->offsets[b+1];++i)
        result(pattern->rows[i],packed)=values(i,v);
    });
    return result;
  }
};

/** @brief Resolve a compact row; generated column rows are sorted exactly. */
inline int column_row(const ColumnTile& tile,int column,int row) {
  const int first=tile.pattern->offsets[column],last=tile.pattern->offsets[column+1];
  if (tile.pattern->full) return first+row;
  const int* begin=tile.pattern->rows.data()+first;
  const int* end=tile.pattern->rows.data()+last;
  const auto found=std::lower_bound(begin,end,row);
  return found!=end && *found==row ? first+(found-begin):-1;
}

/** @brief In-place P_s image on complete resident panels of invariant support.
 * Sigma emits weighted exchange representatives. Their projection has the
 * same bilinear forms against parity-s states, but is itself a parity-s carrier.
 * Each unordered product pair has one owner; only panel pointers are shared.
 */
inline void project_exchange(std::vector<ColumnTile*>& panels,int extent,int parity) {
  profiling::ScopedPhase phase(profiling::Phase::CiProductProjection);
  const int n_beta=(panels.size()-1)*extent+panels.back()->pattern->columns();
  const int workers=coarse_openmp_thread_count(cooperative_openmp_thread_count(),
      static_cast<std::size_t>(n_beta)*n_beta);
  cooperative_openmp_for(n_beta,workers,[&](std::ptrdiff_t b) {
    auto* tile=panels[b/extent];
    if (!tile) return;
    const int column=b%extent,first=tile->pattern->offsets[column];
    for (int i=first;i<tile->pattern->offsets[column+1];++i) {
      const int a=tile->pattern->full ? i-first:tile->pattern->rows[i];
      if (a>b) break;
      auto& other=*panels[a/extent];
      const int j=column_row(other,a%extent,b);
      for (int v=0;v<tile->values.cols();++v) {
        const double value=0.5*(tile->values(i,v)+parity*other.values(j,v));
        tile->values(i,v)=value;
        other.values(j,v)=parity*value;
      }
    }
  });
}

/** @brief Accepted Y/residual panels, retained only when a complete state fits the
 * shared payload budget. Compact support and transform extent are immutable.
 * Larger states use the same jet sweep without retained product data.
 */
struct ProductImages {
  int extent=0;
  std::size_t bytes=0;
  std::size_t budget=workspace::available()/4;
  workspace::Reservation residency;
  std::vector<std::vector<ColumnTile>> states;
  /** Single-root density already contracted for the energy-multiplier term. */
  Eigen::MatrixXd density;
};

/** @brief Gather borrowed columns directly onto exact structural coordinates. */
template<class Column>
ColumnTile gather_columns(int n_alpha,int width,Column&& column,
    std::shared_ptr<const ColumnPattern> pattern) {
  const int cols=pattern->columns();
  ColumnTile result{pattern,Eigen::MatrixXd(pattern->offsets[cols],width)};
  const int workers=coarse_openmp_thread_count(cooperative_openmp_thread_count(),
      static_cast<std::size_t>(n_alpha)*cols*width);
  cooperative_openmp_for(cols*width,workers,[&](std::ptrdiff_t packed) {
    const int b=packed%cols,v=packed/cols;
    const auto source=column(b,v);
    if (pattern->full) {
      result.values.col(v).segment(b*n_alpha,n_alpha)=source;
      return;
    }
    for (int i=pattern->offsets[b];i<pattern->offsets[b+1];++i)
      result.values(i,v)=source[pattern->rows[i]];
  });
  return result;
}

/** @brief Pack structural coordinates, never an amplitude-thresholded sparse view. */
inline ColumnTile pack_columns(Eigen::MatrixXd dense,std::shared_ptr<const ColumnPattern> pattern) {
  profiling::ScopedPhase phase(profiling::Phase::CiProductPack);
  const int cols=pattern->columns(),width=dense.cols()/cols;
  if (pattern->full) {
    dense.resize(pattern->offsets[cols],width);
    return {std::move(pattern),std::move(dense)};
  }
  return gather_columns(dense.rows(),width,[&](int b,int v) {
    return dense.col(v*cols+b);
  },std::move(pattern));
}

/** @brief Pack independent RHS panels without a concatenated dense temporary.
 * Each panel has n_alpha rows and pattern->columns() columns. Only pointer
 * metadata is borrowed; the common gather finishes before these owners die.
 */
inline ColumnTile pack_columns(std::initializer_list<const Eigen::MatrixXd*> panels,
    std::shared_ptr<const ColumnPattern> pattern) {
  profiling::ScopedPhase phase(profiling::Phase::CiProductPack);
  return gather_columns((*panels.begin())->rows(),panels.size(),[&](int b,int v) {
    return panels.begin()[v]->col(b);
  },std::move(pattern));
}

/** @brief Bounded compact panels; eviction consumes partial images linearly. */
class ColumnPanels {
public:
  using Build=std::function<ColumnTile(int)>;
  using Consume=std::function<void(int,ColumnTile&)>;

  ColumnPanels(std::size_t bytes,Build build,Consume consume={})
      : budget_(std::max<std::size_t>(1,bytes)),build_(std::move(build)),consume_(std::move(consume)) {}

  bool contains(int index) const { return index_.contains(index); }

  ColumnTile& get(int index) {
    const auto found=index_.find(index);
    if (found!=index_.end()) {
      panels_.splice(panels_.begin(),panels_,found->second);
      return found->second->tile;
    }
    return insert(index,build_(index));
  }

  /** @brief Admit a newly generated tile; no existing owner may have this key. */
  ColumnTile& insert(int index,ColumnTile tile) {
    const auto bytes=tile.bytes();
    while (!panels_.empty() && (bytes>budget_ || resident_>budget_-bytes)) {
      auto& last=panels_.back();
      resident_-=last.tile.bytes();
      if (consume_) consume_(last.index,last.tile);
      index_.erase(last.index);
      panels_.pop_back();
      ++evictions;
    }
    panels_.push_front({index,std::move(tile)});
    resident_+=bytes;
    peak_bytes=std::max(peak_bytes,resident_);
    index_[index]=panels_.begin();
    return panels_.front().tile;
  }

  /** @brief Consume and free a column whose final physical user has completed. */
  void release(int index) {
    const auto found=index_.find(index);
    if (found==index_.end()) return;
    auto panel=found->second;
    resident_-=panel->tile.bytes();
    if (consume_) consume_(panel->index,panel->tile);
    panels_.erase(panel);
    index_.erase(found);
    ++releases;
  }

  void finish() {
    if (consume_) for (auto& panel:panels_) consume_(panel.index,panel.tile);
    panels_.clear(); index_.clear(); resident_=0;
  }

  std::size_t evictions=0,peak_bytes=0,releases=0;

private:
  struct Entry { int index; ColumnTile tile; };
  std::size_t budget_,resident_=0;
  Build build_;
  Consume consume_;
  std::list<Entry> panels_;
  std::unordered_map<int,std::list<Entry>::iterator> index_;
};

} // namespace xmvb::vb::detail
