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
      : n_orbitals_(n_orbitals), alpha_(prefixes(alpha, n_orbitals)),
        beta_(prefixes(beta, n_orbitals)) {}

  bool contains(int alpha, int beta, int depth = 0) const noexcept {
    const auto* a = alpha_.data() + static_cast<std::size_t>(alpha) * n_orbitals_;
    const auto* b = beta_.data() + static_cast<std::size_t>(beta) * n_orbitals_;
    for (int k = 0; k < n_orbitals_; ++k)
      if (static_cast<int>(a[k]) + b[k] + depth < k + 1) return false;
    return true;
  }

private:
  static std::vector<std::uint8_t> prefixes(
      const std::vector<std::uint64_t>& masks, int n_orbitals) {
    std::vector<std::uint8_t> result(masks.size() * n_orbitals);
    for (std::size_t string = 0; string < masks.size(); ++string) {
      int count = 0;
      for (int k = 0; k < n_orbitals; ++k) {
        count += (masks[string] >> k) & 1ULL;
        result[string * n_orbitals + k] = static_cast<std::uint8_t>(count);
      }
    }
    return result;
  }
  int n_orbitals_;
  std::vector<std::uint8_t> alpha_, beta_;
};

}  // namespace xmvb::mcscf
