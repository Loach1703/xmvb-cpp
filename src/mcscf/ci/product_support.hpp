#pragma once

#include <algorithm>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

#include "mcscf/ci/covalent_support.hpp"

namespace xmvb::mcscf {

/** @brief Downward closure of selected alpha/beta products and its generator halo.
 * A product is reachable from seed (A,B) when each spin's prefix counts
 * dominate that seed. One generator may reduce the counts by one in ONE
 * spin sector. A tree of prefix lower bounds accelerates the union query;
 * storage is O((U_alpha+U_beta+N_seed)*n_orbitals), never a product mask.
 */
class ProductSupport {
public:
  /** @brief Complete singly occupied seed family (the covalent closed form). */
  ProductSupport(int n_orbitals, const std::vector<std::uint64_t>& alpha,
                 const std::vector<std::uint64_t>& beta)
      : chunks_((n_orbitals + 7) / 8), n_beta_(beta.size()),
        covalent_(std::in_place, n_orbitals, alpha, beta), complete_covalent_(true) {}

  /** @brief Seeds are flattened alpha*n_beta+beta indices, not CI amplitudes. */
  ProductSupport(int n_orbitals, const std::vector<std::uint64_t>& alpha,
                 const std::vector<std::uint64_t>& beta,
                 const std::vector<int>& products)
      : chunks_((n_orbitals + 7) / 8), n_beta_(beta.size()) {
    const auto full = (std::uint64_t{1} << n_orbitals) - 1;
    complete_covalent_ = products.size() == alpha.size() &&
        products.size() == beta.size();
    seeds_.reserve(products.size());
    for (const int product : products) {
      if (product < 0 || static_cast<std::size_t>(product) >= alpha.size()*beta.size())
        throw std::invalid_argument("product support seed is outside the spin spaces");
      seeds_.push_back({product / n_beta_, product % n_beta_});
      const auto a=alpha[product/n_beta_], b=beta[product%n_beta_];
      const int na=__builtin_popcountll(a), nb=__builtin_popcountll(b);
      full_ = full_ || (a == (full ^ ((std::uint64_t{1} << (n_orbitals-na))-1)) &&
                       b == (full ^ ((std::uint64_t{1} << (n_orbitals-nb))-1)));
      complete_covalent_ = complete_covalent_ &&
          (alpha[product / n_beta_] & beta[product % n_beta_]) == 0 &&
          (alpha[product / n_beta_] | beta[product % n_beta_]) == full;
    }
    std::sort(seeds_.begin(), seeds_.end(), [](const Seed& a, const Seed& b) {
      return a.alpha != b.alpha ? a.alpha < b.alpha : a.beta < b.beta;
    });
    seeds_.erase(std::unique(seeds_.begin(), seeds_.end(), [](const Seed& a, const Seed& b) {
      return a.alpha == b.alpha && a.beta == b.beta;
    }), seeds_.end());
    complete_covalent_ = complete_covalent_ && seeds_.size() == alpha.size();
    if (complete_covalent_) covalent_.emplace(n_orbitals,alpha,beta);
    else if (!full_ && !seeds_.empty()) {
      alpha_=prefixes(alpha,n_orbitals);
      beta_=prefixes(beta,n_orbitals);
      build(0,seeds_.size());
    }
    if (full_ || complete_covalent_) std::vector<Seed>().swap(seeds_);
  }

  bool full() const noexcept { return full_; }
  bool contains(int alpha, int beta, int depth = 0) const noexcept {
    if (full_) return true;
    if (complete_covalent_) return covalent_->contains(alpha, beta, depth);
    return !nodes_.empty() && contains_node(0,
        alpha_.data() + static_cast<std::size_t>(alpha)*chunks_,
        beta_.data() + static_cast<std::size_t>(beta)*chunks_, depth);
  }

  std::size_t dynamic_bytes() const noexcept {
    return (alpha_.capacity()+beta_.capacity()+bounds_.capacity())*sizeof(std::uint64_t)+
        seeds_.capacity()*sizeof(Seed)+nodes_.capacity()*sizeof(Node)+
        (covalent_ ? covalent_->dynamic_bytes() : 0);
  }

private:
  struct Seed { int alpha, beta; };
  struct Node { int first, last, left = -1, right = -1; };
  static constexpr std::uint64_t high_ = 0x8080808080808080ULL;
  static constexpr std::uint64_t one_ = 0x0101010101010101ULL;

  static std::vector<std::uint64_t> prefixes(
      const std::vector<std::uint64_t>& masks, int n) {
    const int chunks = (n+7)/8;
    std::vector<std::uint64_t> result(masks.size()*chunks, 0);
    for (std::size_t i=0; i<masks.size(); ++i) {
      int count=0;
      for (int k=0; k<n; ++k) {
        count += (masks[i] >> k)&1ULL;
        result[i*chunks+k/8] |= std::uint64_t(count) << (8*(k%8));
      }
    }
    return result;
  }
  bool dominates(const std::uint64_t* value, const std::uint64_t* bound,
                 int depth) const noexcept {
    for (int k=0; k<chunks_; ++k)
      if (((((value[k] + depth*one_) | high_) - bound[k]) & high_) != high_)
        return false;
    return true;
  }
  bool matches(const std::uint64_t* a, const std::uint64_t* b,
               const std::uint64_t* sa, const std::uint64_t* sb,
               int depth) const noexcept {
    const bool a0 = dominates(a, sa, 0);
    if (a0 && dominates(b, sb, depth)) return true;
    return depth != 0 && !a0 && dominates(a, sa, 1) && dominates(b, sb, 0);
  }
  bool contains_node(int index, const std::uint64_t* a,
                     const std::uint64_t* b, int depth) const noexcept {
    const auto* bound = bounds_.data() + static_cast<std::size_t>(index)*2*chunks_;
    if (!matches(a,b,bound,bound+chunks_,depth)) return false;
    const Node& node = nodes_[index];
    if (node.left >= 0)
      return contains_node(node.left,a,b,depth) || contains_node(node.right,a,b,depth);
    for (int i=node.first; i<node.last; ++i)
      if (matches(a,b,alpha_.data()+static_cast<std::size_t>(seeds_[i].alpha)*chunks_,
          beta_.data()+static_cast<std::size_t>(seeds_[i].beta)*chunks_,depth)) return true;
    return false;
  }
  int build(int first, int last) {
    const int index=nodes_.size();
    nodes_.push_back({first,last});
    const std::size_t offset=bounds_.size();
    bounds_.resize(offset+2*chunks_, 0);
    int split_lane=0, largest_span=-1;
    for (int lane=0; lane<16*chunks_; ++lane) {
      int low=255, high=0;
      for (int i=first; i<last; ++i) {
        const int value=prefix_lane(seeds_[i],lane);
        low=std::min(low,value); high=std::max(high,value);
      }
      bounds_[offset+lane/8] |= std::uint64_t(low) << (8*(lane%8));
      if (high-low > largest_span) { largest_span=high-low; split_lane=lane; }
    }
    if (last-first > 8) {
      const int middle=(first+last)/2;
      std::nth_element(seeds_.begin()+first,seeds_.begin()+middle,seeds_.begin()+last,
          [&](const Seed& a,const Seed& b) {return prefix_lane(a,split_lane)<prefix_lane(b,split_lane);});
      const int left=build(first,middle), right=build(middle,last);
      nodes_[index].left=left; nodes_[index].right=right;
    }
    return index;
  }
  int prefix_lane(const Seed& seed, int lane) const noexcept {
    const bool alpha=lane<8*chunks_;
    const int k=alpha?lane:lane-8*chunks_;
    const auto& counts=alpha?alpha_:beta_;
    const int string=alpha?seed.alpha:seed.beta;
    return (counts[static_cast<std::size_t>(string)*chunks_+k/8] >> (8*(k%8)))&255ULL;
  }
  int chunks_, n_beta_;
  std::optional<CovalentSupport> covalent_;
  bool complete_covalent_=false;
  bool full_=false;
  std::vector<std::uint64_t> alpha_, beta_, bounds_;
  std::vector<Seed> seeds_;
  std::vector<Node> nodes_;
};

} // namespace xmvb::mcscf
