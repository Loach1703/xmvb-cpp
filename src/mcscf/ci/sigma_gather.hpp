#pragma once

#include <algorithm>
#include <span>

namespace xmvb::mcscf::detail {

/** @brief One tile-local target with its composed grouped-image address.
 * An RHS advances by stride image rows. All temporal responses share this
 * symbolic address; no channel values or excitation-edge arrays are cached.
 */
struct SigmaGatherRow { int alpha,image,stride; };

/** @brief Contract one beta edge over existing target rows only.
 * Targets are sorted alpha addresses with the layout permutation composed once.
 * Sparse destination rows are intersected exactly, without a product mask.
 * The caller owns a disjoint destination interval, so no atomics are needed.
 */
inline void gather_beta_rows(const double* values,int pairs,int rhs,double sign,
    std::span<const SigmaGatherRow> targets,int first,int last,
    double* destination,const int* rows,int n_rows) {
  const auto before=[](const SigmaGatherRow& row,int alpha) { return row.alpha<alpha; };
  auto target=std::lower_bound(targets.begin(),targets.end(),first,before);
  const auto end=std::lower_bound(target,targets.end(),last,before);
  if (target==end) return;
  if (!rows) {
    for (;target!=end;++target)
      destination[target->alpha]+=sign*values[
          static_cast<std::size_t>(target->image+rhs*target->stride)*pairs];
    return;
  }
  int output=std::lower_bound(rows,rows+n_rows,target->alpha)-rows;
  for (;target!=end;++target) {
    const int alpha=target->alpha;
    while (output<n_rows && rows[output]<alpha) ++output;
    if (output==n_rows) break;
    if (rows[output]==alpha)
      destination[output]+=sign*values[
          static_cast<std::size_t>(target->image+rhs*target->stride)*pairs];
  }
}

} // namespace xmvb::mcscf::detail
