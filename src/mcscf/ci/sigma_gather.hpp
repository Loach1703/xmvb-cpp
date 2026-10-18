#pragma once

#include <algorithm>
#include <span>

#include "mcscf/ci/channel_product.hpp"

namespace xmvb::mcscf::detail {

/** @brief Contract one beta edge over existing target rows only.
 * Targets are sorted halo addresses. Their index entries contain original
 * tile rows; the shared layout maps each RHS to its grouped image.
 * Sparse destination rows are intersected exactly, without a product mask.
 * The caller owns a disjoint destination interval, so no atomics are needed.
 */
inline void gather_beta_rows(const double* values,int pairs,int rhs_offset,double sign,
    std::span<const int> targets,const int* index,std::span<const int> halo,
    const ChannelLayout& layout,int first,int last,double* destination,const int* rows,int n_rows) {
  auto target=std::lower_bound(targets.begin(),targets.end(),first);
  const auto end=std::lower_bound(target,targets.end(),last);
  if (target==end) return;
  if (!rows) {
    for (;target!=end;++target)
      destination[halo[*target]]+=sign*values[
          static_cast<std::size_t>(layout.image(index[*target]+rhs_offset))*pairs];
    return;
  }
  int output=std::lower_bound(rows,rows+n_rows,halo[*target])-rows;
  for (;target!=end;++target) {
    const int alpha=halo[*target];
    while (output<n_rows && rows[output]<alpha) ++output;
    if (output==n_rows) break;
    if (rows[output]==alpha)
      destination[output]+=sign*values[
          static_cast<std::size_t>(layout.image(index[*target]+rhs_offset))*pairs];
  }
}

} // namespace xmvb::mcscf::detail
