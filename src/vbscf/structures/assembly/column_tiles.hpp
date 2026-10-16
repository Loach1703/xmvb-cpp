#pragma once

#include <algorithm>
#include <functional>
#include <list>
#include <memory>
#include <unordered_map>

#include <Eigen/Core>

#include "core/openmp.hpp"

namespace xmvb::vb::detail {

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

/** @brief Pack structural coordinates, never an amplitude-thresholded sparse view. */
inline ColumnTile pack_columns(Eigen::MatrixXd dense,std::shared_ptr<const ColumnPattern> pattern) {
  const int cols=pattern->columns(),width=dense.cols()/cols;
  if (pattern->full) {
    dense.resize(pattern->offsets[cols],width);
    return {std::move(pattern),std::move(dense)};
  }
  ColumnTile result{pattern,Eigen::MatrixXd(pattern->offsets[cols],width)};
  const int workers=coarse_openmp_thread_count(cooperative_openmp_thread_count(),dense.size());
  cooperative_openmp_for(cols*width,workers,[&](std::ptrdiff_t packed) {
    const int b=packed%cols,v=packed/cols;
    for (int i=pattern->offsets[b];i<pattern->offsets[b+1];++i)
      result.values(i,v)=dense(pattern->rows[i],packed);
  });
  return result;
}

/** @brief Bounded compact panels; eviction consumes partial images linearly. */
class ColumnPanels {
public:
  using Build=std::function<ColumnTile(int)>;
  using Consume=std::function<void(int,ColumnTile&)>;

  ColumnPanels(std::size_t bytes,Build build,Consume consume={})
      : budget_(std::max<std::size_t>(1,bytes)),build_(std::move(build)),consume_(std::move(consume)) {}

  ColumnTile& get(int index) {
    const auto found=index_.find(index);
    if (found!=index_.end()) {
      panels_.splice(panels_.begin(),panels_,found->second);
      return found->second->tile;
    }
    ColumnTile tile=build_(index);
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

  void finish() {
    if (consume_) for (auto& panel:panels_) consume_(panel.index,panel.tile);
    panels_.clear(); index_.clear(); resident_=0;
  }

  std::size_t evictions=0,peak_bytes=0;

private:
  struct Entry { int index; ColumnTile tile; };
  std::size_t budget_,resident_=0;
  Build build_;
  Consume consume_;
  std::list<Entry> panels_;
  std::unordered_map<int,std::list<Entry>::iterator> index_;
};

} // namespace xmvb::vb::detail
