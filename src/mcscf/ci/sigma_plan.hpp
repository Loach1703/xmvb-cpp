#pragma once

#include <cstddef>
#include <optional>
#include <vector>

#include "mcscf/ci/links.hpp"
#include "mcscf/ci/product_support.hpp"

namespace xmvb::mcscf::detail {

/** @brief Numeric halo, channel and target-index budget; immutable spin links are shared. */
struct SigmaTileBudget {
  int width,pairs,slabs;
  std::size_t bytes;
  int packed_slabs=0;
};

struct SigmaTargetGroup { int beta,begin,end,column=0; };

/** @brief Stream each target once, admitting its actual shared neighborhoods.
 * One pending target crosses a budget boundary without rewinding the support
 * cursor. Storage is single-spin metadata plus the current tile, never U^2.
 */
class SigmaTilePlan {
public:
  SigmaTilePlan(const SpinConnections& alpha,const SpinConnections& beta,
      const ProductSupport* support,int parity,SigmaTileBudget budget,bool release);
  bool next();
  std::size_t payload_bytes() const;

  std::vector<SigmaTargetGroup> groups;
  std::vector<int> active,group_index,halo,columns,halo_index,column_index,released;
  std::vector<unsigned char> active_support;
  std::vector<int> channel_groups;

private:
  struct Node { int alpha,beta; bool supported; };
  bool pull(Node& node);
  long double payload(std::size_t rows,std::size_t h,std::size_t b,std::size_t g,
      std::size_t channels) const;
  const SpinConnections& alpha_;
  const SpinConnections& beta_;
  const ProductSupport* support_;
  int parity_;
  SigmaTileBudget budget_;
  std::vector<int> beta_targets_,seen_alpha_,new_alpha_,new_beta_,last_use_,release_order_;
  std::vector<int> column_rows_;
  std::vector<unsigned char> alpha_seen_;
  std::optional<ProductSupport::Cursor> cursor_;
  std::optional<Node> pending_;
  std::size_t next_beta_=0,next_release_=0;
  std::size_t channels_=0;
  int next_alpha_=0;
};

} // namespace xmvb::mcscf::detail
