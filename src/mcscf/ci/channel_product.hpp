#pragma once

#include <algorithm>
#include <span>
#include <vector>
#include <Eigen/Core>

namespace xmvb::mcscf::detail {

using ChannelRows=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;

/** @brief Exact source-channel restriction; every output channel is retained.
 * Group zero is unrestricted. A nonzero group 1+first*n+last permits packed
 * pairs (p,q) with p<=first<=last<q. Small immutable kernel panels are shared
 * by workers. Packing borrows disjoint slices of the consumed input halo,
 * so no large per-thread allocation or additional channel slab is needed.
 */
class ChannelProduct {
  struct Block { std::vector<int> columns; Eigen::MatrixXd kernel; };
  const Eigen::MatrixXd& kernel_;
  std::vector<Block> blocks_;
public:
  ChannelProduct(const Eigen::MatrixXd& kernel,int orbitals,bool restricted):kernel_(kernel) {
    if (!restricted) return;
    blocks_.resize(1+orbitals*orbitals);
    for (int first=0;first<orbitals;++first) for (int last=first;last<orbitals-1;++last) {
      auto& block=blocks_[1+first*orbitals+last];
      for (int q=last+1;q<orbitals;++q) for (int p=0;p<=first;++p)
        block.columns.push_back(q*(q+1)/2+p);
      block.kernel.resize(kernel.rows(),block.columns.size());
      for (int k=0;k<static_cast<int>(block.columns.size());++k)
        block.kernel.col(k)=kernel.col(block.columns[k]);
    }
  }

  /** @brief Thread-owned row interval. Call after every reader consumed halo.
   * Source/destination may alias: all channels of a row are packed before
   * that row is overwritten. Dense in-place work also uses this scratch.
   */
  void apply(Eigen::Ref<const ChannelRows> source,Eigen::Ref<ChannelRows> destination,
      int first,int last,std::span<const int> groups,double* scratch,std::size_t size) const {
    if (first==last) return;
    const int pairs=kernel_.rows();
    if (blocks_.empty() || groups.empty() || size<static_cast<std::size_t>(2*pairs)) {
      if (source.data()!=destination.data())
        destination.middleRows(first,last-first).noalias()=
            source.middleRows(first,last-first)*kernel_.transpose();
      else for (int row=first;row<last;++row)
        destination.row(row)=(source.row(row)*kernel_.transpose()).eval();
      return;
    }
    std::vector<int> offsets(blocks_.size()+1,0),order(last-first);
    for (int row=first;row<last;++row) ++offsets[groups[row%groups.size()]+1];
    for (std::size_t k=1;k<offsets.size();++k) offsets[k]+=offsets[k-1];
    auto position=offsets;
    for (int row=first;row<last;++row) order[position[groups[row%groups.size()]]++]=row;
    for (int group=0;group<static_cast<int>(blocks_.size());++group) {
      const auto& block=blocks_[group];
      const int inputs=group ? block.columns.size():pairs;
      const int extent=std::min<std::size_t>(last-first,size/(inputs+pairs));
      for (int begin=offsets[group];begin<offsets[group+1];begin+=extent) {
        const int rows=std::min(extent,offsets[group+1]-begin);
        Eigen::Map<ChannelRows> packed(scratch,rows,inputs);
        Eigen::Map<ChannelRows> image(scratch+static_cast<std::size_t>(rows)*inputs,rows,pairs);
        for (int i=0;i<rows;++i) {
          const int row=order[begin+i];
          if (!group) packed.row(i)=source.row(row);
          else for (int k=0;k<inputs;++k) packed(i,k)=source(row,block.columns[k]);
        }
        if (group) image.noalias()=packed*block.kernel.transpose();
        else image.noalias()=packed*kernel_.transpose();
        for (int i=0;i<rows;++i) destination.row(order[begin+i])=image.row(i);
      }
    }
  }
};

} // namespace xmvb::mcscf::detail
