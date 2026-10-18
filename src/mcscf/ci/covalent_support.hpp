#pragma once

#include <cstdint>
#include <bit>
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace xmvb::mcscf {

/** @brief Exact support of an upper-triangular transform of singly occupied sites.
 * Prefix electron counts must be at least the prefix length. A single
 * one-body excitation can lower a prefix by one: depth=1 is its exact
 * structural closure. Metadata is linear in the number of spin strings;
 * no product-space mask is retained here.
 */
class CovalentSupport {
  enum Edge : std::uint8_t { Empty=1, Occupied=2, Complete=4 };
public:
  CovalentSupport(int n_orbitals, const std::vector<std::uint64_t>& alpha,
                  const std::vector<std::uint64_t>& beta)
      : chunks_((n_orbitals + 7) / 8), alpha_(prefixes(alpha, n_orbitals)),
        beta_(prefixes(beta, n_orbitals)), lengths_(chunks_, 0),
        n_(n_orbitals), electrons_(std::popcount(alpha.front())), stride_(electrons_+1),
        choose_((n_+1)*stride_,0), branches_(beta.size()*(n_+1)*stride_,0),
        alpha_order_(alpha.size()) {
    for (int k = 0; k < n_orbitals; ++k)
      lengths_[k / 8] |= static_cast<std::uint64_t>(k + 1) << (8 * (k % 8));
    // C(k,e) also gives each complete subtree's colex rank interval.
    for (int k=0;k<=n_;++k) {
      choose_[k*stride_]=1;
      for (int e=1;e<=std::min(k,electrons_);++e)
        choose_[k*stride_+e]=choose_[(k-1)*stride_+e-1]+choose_[(k-1)*stride_+e];
    }
    if (choose_[n_*stride_+electrons_]!=alpha.size())
      throw std::invalid_argument("covalent column plan requires complete alpha strings");
    for (int i=0;i<static_cast<int>(alpha.size());++i) {
      std::size_t rank=0;
      int e=0;
      for (auto mask=alpha[i];mask;mask&=mask-1)
        rank+=choose_[std::countr_zero(mask)*stride_+(++e)];
      alpha_order_[rank]=i;
    }
    ordered_=std::is_sorted(alpha_order_.begin(),alpha_order_.end());
    // A node (k,e) represents valid alpha occupations on the first k sites.
    // Its two edges and complete-subtree flag depend on beta, never on a CI
    // RHS or orbital coefficients. Retain only three bits, not product rows.
    std::vector<std::uint64_t> counts(choose_.size());
    for (std::size_t b=0;b<beta.size();++b) {
      std::fill(counts.begin(),counts.end(),0);
      counts[0]=1;
      auto* nodes=branches_.data()+b*choose_.size();
      nodes[0]=Complete;
      int beta_count=0;
      for (int k=1;k<=n_;++k) {
        beta_count+=(beta[b]>>(k-1))&1;
        for (int e=0;e<=std::min(k,electrons_);++e) {
          if (e+beta_count<k) continue;
          const auto zero=counts[(k-1)*stride_+e];
          const auto one=e ? counts[(k-1)*stride_+e-1]:0;
          counts[k*stride_+e]=zero+one;
          nodes[k*stride_+e]=(zero ? Empty:0)|(one ? Occupied:0)|
              (zero+one==choose_[k*stride_+e] ? Complete:0);
        }
      }
    }
  }

  /** @brief Emit sorted original-alpha indices from a reusable prefix DAG.
   * Entire legal subtrees are colex intervals. Scratch consists only of the
   * bounded traversal stack and the caller's current tile row list. Persistent
   * storage is O(U_beta * n_orbitals * (n_alpha_electrons+1)), never U_alpha*U_beta.
   */
  void append_column_rows(int beta,std::vector<int>& rows) const {
    const auto* nodes=branches_.data()+static_cast<std::size_t>(beta)*choose_.size();
    // Caller indices may be permuted. A current-column bitset restores their
    // order in linear work, without sorting rows or retaining product masks.
    std::vector<std::uint64_t> indices(ordered_ ? 0:(alpha_order_.size()+63)/64,0);
    const auto visit=[&](auto&& self,int k,int e,std::size_t rank)->void {
      const auto edges=nodes[k*stride_+e];
      if (edges&Complete) {
        const auto count=choose_[k*stride_+e];
        if (ordered_) for (std::size_t i=0;i<count;++i) rows.push_back(rank+i);
        else for (std::size_t i=0;i<count;++i) {
          const auto index=alpha_order_[rank+i];
          indices[index/64]|=std::uint64_t{1}<<(index%64);
        }
        return;
      }
      if (edges&Empty) self(self,k-1,e,rank);
      if (edges&Occupied) self(self,k-1,e-1,rank+choose_[(k-1)*stride_+e]);
    };
    visit(visit,n_,electrons_,0);
    for (std::size_t w=0;w<indices.size();++w)
      for (auto bits=indices[w];bits;bits&=bits-1)
        rows.push_back(64*w+std::countr_zero(bits));
  }

  bool contains(int alpha, int beta, int depth = 0) const noexcept {
    const auto* a = alpha_.data() + static_cast<std::size_t>(alpha) * chunks_;
    const auto* b = beta_.data() + static_cast<std::size_t>(beta) * chunks_;
    constexpr std::uint64_t high = 0x8080808080808080ULL;
    constexpr std::uint64_t one = 0x0101010101010101ULL;
    // Eight exact byte-lane comparisons per word. At half filling with
    // <=63 orbitals and depth 0/1, sums stay <128, so lanes cannot carry.
    // The high-bit guard prevents subtraction borrows between byte lanes.
    for (int k = 0; k < chunks_; ++k) {
      const auto sum = a[k] + b[k] + static_cast<std::uint64_t>(depth) * one;
      if ((((sum | high) - lengths_[k]) & high) != high) return false;
    }
    return true;
  }

  /** @brief First/last deficient prefix; {-1,-1} means depth-zero support.
   * A frontier node can be reached only by a raising replacement crossing
   * every deficient prefix. Byte-lane arithmetic shares contains()'s proof.
   */
  std::pair<int,int> deficient_prefixes(int alpha,int beta) const noexcept {
    constexpr std::uint64_t high=0x8080808080808080ULL;
    const auto* a=alpha_.data()+static_cast<std::size_t>(alpha)*chunks_;
    const auto* b=beta_.data()+static_cast<std::size_t>(beta)*chunks_;
    int first=-1,last=-1;
    for (int k=0;k<chunks_;++k) {
      const auto bad=~(((a[k]+b[k])|high)-lengths_[k])&high;
      if (bad) {
        if (first<0) first=8*k+std::countr_zero(bad)/8;
        last=8*k+(63-std::countl_zero(bad))/8;
      }
    }
    return {first,last};
  }

  std::size_t dynamic_bytes() const noexcept {
    return (alpha_.capacity()+beta_.capacity()+lengths_.capacity()+choose_.capacity())*
        sizeof(std::uint64_t)+branches_.capacity()*sizeof(std::uint8_t)+
        alpha_order_.capacity()*sizeof(int);
  }

private:
  static std::vector<std::uint64_t> prefixes(
      const std::vector<std::uint64_t>& masks, int n_orbitals) {
    const int chunks = (n_orbitals + 7) / 8;
    std::vector<std::uint64_t> result(masks.size() * chunks, 0);
    for (std::size_t string = 0; string < masks.size(); ++string) {
      int count = 0;
      for (int k = 0; k < n_orbitals; ++k) {
        count += (masks[string] >> k) & 1ULL;
        result[string * chunks + k / 8] |=
            static_cast<std::uint64_t>(count) << (8 * (k % 8));
      }
    }
    return result;
  }
  int chunks_;
  std::vector<std::uint64_t> alpha_, beta_, lengths_;
  int n_,electrons_,stride_;
  std::vector<std::uint64_t> choose_;
  std::vector<std::uint8_t> branches_;
  std::vector<int> alpha_order_;
  bool ordered_=false;
};

}  // namespace xmvb::mcscf
