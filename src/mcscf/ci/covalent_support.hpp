#pragma once

#include <cstdint>
#include <vector>

namespace xmvb::mcscf {

/** @brief Exact support of an upper-triangular transform of singly occupied sites.
 * Prefix electron counts must be at least the prefix length. A single
 * one-body excitation can lower a prefix by one: depth=1 is its exact
 * structural closure. Metadata is linear in the number of spin strings;
 * no product-space mask is retained here.
 */
class CovalentSupport {
public:
  CovalentSupport(int n_orbitals, const std::vector<std::uint64_t>& alpha,
                  const std::vector<std::uint64_t>& beta)
      : chunks_((n_orbitals + 7) / 8), alpha_(prefixes(alpha, n_orbitals)),
        beta_(prefixes(beta, n_orbitals)), lengths_(chunks_, 0) {
    for (int k = 0; k < n_orbitals; ++k)
      lengths_[k / 8] |= static_cast<std::uint64_t>(k + 1) << (8 * (k % 8));
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
};

}  // namespace xmvb::mcscf
