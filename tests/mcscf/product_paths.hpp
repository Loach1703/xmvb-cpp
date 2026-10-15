#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace xmvb::test {

/** @brief Exact prefix-path prototype; no product carrier or product mask.
 * A seed is an original touched alpha/beta spin-string pair, not a VB
 * structure. Multiple structures contribute to the same seed before this
 * traversal. The callback receives each reachable product exactly once.
 * Kept in tests until the production contraction scheduler consumes paths.
 */
class ProductPaths {
public:
  using Mask = std::uint64_t;
  using Seed = std::pair<Mask, Mask>;
  struct Counts {
    std::size_t prefixes = 0;
    std::size_t words = 0;
    std::size_t products = 0;
    std::size_t support = 0;
  };

  ProductPaths(int n, int n_alpha, int n_beta, std::vector<Seed> seeds)
      : n_(n), na_(n_alpha), nb_(n_beta) {
    if (n < 1 || n > 63 || n_alpha < 0 || n_alpha > n ||
        n_beta < 0 || n_beta > n)
      throw std::invalid_argument("invalid prefix-path spin sector");
    const Mask full = (Mask{1} << n) - 1;
    for (const auto& [a, b] : seeds)
      if ((a & ~full) || (b & ~full) || __builtin_popcountll(a) != na_ ||
          __builtin_popcountll(b) != nb_)
        throw std::invalid_argument("prefix-path seed is outside the spin sector");
    std::sort(seeds.begin(), seeds.end());
    seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());
    size_ = seeds.size();
    words_ = (size_ + 63) / 64;
    bounds_.assign(2 * static_cast<std::size_t>(n_) * (n_ + 1) * words_, 0);
    for (std::size_t seed = 0; seed < size_; ++seed) {
      int a = 0, b = 0;
      for (int k = 0; k < n_; ++k) {
        a += (seeds[seed].first >> k) & 1;
        b += (seeds[seed].second >> k) & 1;
        bounds_[offset(0, k, a) + seed / 64] |= Mask{1} << (seed % 64);
        bounds_[offset(1, k, b) + seed / 64] |= Mask{1} << (seed % 64);
      }
    }
    // L_sigma(k,c) contains exactly seeds with p_k(seed_sigma) <= c.
    for (int spin = 0; spin < 2; ++spin)
      for (int k = 0; k < n_; ++k)
        for (int count = 1; count <= n_; ++count)
          for (std::size_t w = 0; w < words_; ++w)
            bounds_[offset(spin, k, count) + w] |=
                bounds_[offset(spin, k, count - 1) + w];
  }

  /** @brief Streams S (depth=0) or its one-generator closure I (depth=1).
   * The third callback argument marks membership in S, even when visiting I.
   * The three seed sets keep the same seed in both spin sectors; the depth
   * allowance cannot switch spin sectors halfway along a path.
   */
  template<class Emit>
  Counts visit(int depth, Emit&& emit) const {
    if (depth != 0 && depth != 1)
      throw std::invalid_argument("prefix-path depth must be zero or one");
    Counts counts;
    if (size_ == 0) return counts;
    const int modes = depth == 0 ? 1 : 3;
    const std::size_t stride = modes * words_;
    std::vector<Mask> live((n_ + 1) * stride, 0);
    for (int mode = 0; mode < modes; ++mode) {
      std::fill_n(live.data() + mode * words_, words_, ~Mask{0});
      if (size_ % 64)
        live[mode * words_ + words_ - 1] = (Mask{1} << (size_ % 64)) - 1;
    }
    const auto descend = [&](const auto& self, int k, int a, int b,
                             Mask alpha, Mask beta, bool in_support) -> void {
      ++counts.prefixes;
      if (k == n_) {
        ++counts.products;
        counts.support += in_support;
        emit(alpha, beta, in_support);
        return;
      }
      const Mask* parent = live.data() + k * stride;
      Mask* child = live.data() + (k + 1) * stride;
      for (int da = 0; da <= 1; ++da) {
        const int next_a = a + da;
        if (next_a > na_ || next_a + n_ - k - 1 < na_) continue;
        for (int db = 0; db <= 1; ++db) {
          const int next_b = b + db;
          if (next_b > nb_ || next_b + n_ - k - 1 < nb_) continue;
          bool any = false, zero = false;
          for (int mode = 0; mode < modes; ++mode) {
            const auto a_offset = offset(0, k, std::min(n_, next_a + (mode == 1)));
            const auto b_offset = offset(1, k, std::min(n_, next_b + (mode == 2)));
            Mask present = 0;
            for (std::size_t w = 0; w < words_; ++w) {
              const Mask value = parent[mode * words_ + w] &
                  bounds_[a_offset + w] & bounds_[b_offset + w];
              child[mode * words_ + w] = value;
              present |= value;
            }
            counts.words += words_;
            any = any || present != 0;
            if (mode == 0) zero = present != 0;
          }
          if (any)
            self(self, k + 1, next_a, next_b,
                 alpha | (static_cast<Mask>(da) << k),
                 beta | (static_cast<Mask>(db) << k), zero);
        }
      }
    };
    descend(descend, 0, 0, 0, 0, 0, true);
    return counts;
  }

  std::size_t metadata_bytes() const noexcept { return bounds_.capacity() * sizeof(Mask); }
  std::size_t workspace_bytes(int depth) const noexcept {
    return (n_ + 1) * (depth == 0 ? 1 : 3) * words_ * sizeof(Mask);
  }

private:
  std::size_t offset(int spin, int k, int count) const noexcept {
    return ((static_cast<std::size_t>(spin) * n_ + k) * (n_ + 1) + count) * words_;
  }
  int n_, na_, nb_;
  std::size_t size_ = 0, words_ = 0;
  std::vector<Mask> bounds_;
};

} // namespace xmvb::test
