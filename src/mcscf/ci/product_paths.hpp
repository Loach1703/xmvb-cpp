#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace xmvb::mcscf {

/** @brief Direct, resumable enumeration of selected-product support.
 *
 * Seeds are distinct original alpha/beta product pairs, after contributions
 * from all selected structures to an original product have been combined.
 * A depth-one path spends its allowance in one fixed spin sector. The cursor
 * owns only traversal stacks and seed-index bitsets; it never stores the
 * emitted products or constructs a complete spin-string universe.
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

  /** @param complete_covalent Certifies the complete singly occupied seed family.
   * This certificate requires n_alpha+n_beta=n and may be supplied without
   * materializing that family in seeds.
   */
  ProductPaths(int n, int n_alpha, int n_beta, std::vector<Seed> seeds,
               bool complete_covalent = false)
      : n_(n), na_(n_alpha), nb_(n_beta), covalent_(complete_covalent) {
    if (n < 1 || n > 63 || na_ < 0 || na_ > n || nb_ < 0 || nb_ > n ||
        (covalent_ && na_ + nb_ != n_))
      throw std::invalid_argument("invalid prefix-path spin sector");
    const Mask full = (Mask{1} << n_) - 1;
    for (const auto& [a, b] : seeds)
      if ((a & ~full) || (b & ~full) || __builtin_popcountll(a) != na_ ||
          __builtin_popcountll(b) != nb_)
        throw std::invalid_argument("prefix-path seed is outside the spin sector");
    if (covalent_) return;

    // The maximal pair dominates every product in the fixed spin sector.
    const Mask highest_a = full ^ ((Mask{1} << (n_ - na_)) - 1);
    const Mask highest_b = full ^ ((Mask{1} << (n_ - nb_)) - 1);
    if (std::find(seeds.begin(), seeds.end(), Seed{highest_a, highest_b}) !=
        seeds.end()) {
      full_ = true;
      return;
    }
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
    // L_spin(k,c) contains precisely the seeds whose prefix count is <= c.
    for (int spin = 0; spin < 2; ++spin)
      for (int k = 0; k < n_; ++k)
        for (int count = 1; count <= n_; ++count)
          for (std::size_t w = 0; w < words_; ++w)
            bounds_[offset(spin, k, count) + w] |=
                bounds_[offset(spin, k, count - 1) + w];
  }

  /** @brief A pull cursor; the ProductPaths object must outlive it.
   * Each call emits at most one product. Beta masks are ordered by a DFS of
   * their orbital bits, and all alpha masks for one beta are emitted before
   * advancing to the next beta.
   */
  class Cursor {
  public:
    Cursor(const ProductPaths& paths, int depth, bool fixed, Mask beta)
        : paths_(&paths), depth_(depth), fixed_(fixed),
          beta_modes_(depth == 0 ? 1 : 2), alpha_modes_(depth == 0 ? 1 : 3),
          beta_live_((paths.n_ + 1) * beta_modes_ * paths.words_, 0),
          alpha_live_((paths.n_ + 1) * alpha_modes_ * paths.words_, 0) {
      if (depth != 0 && depth != 1)
        throw std::invalid_argument("prefix-path depth must be zero or one");
      beta_stack_.reserve(paths.n_ + 1);
      alpha_stack_.reserve(paths.n_ + 1);
      if (fixed && ((beta >> paths.n_) || __builtin_popcountll(beta) != paths.nb_))
        throw std::invalid_argument("prefix-path beta is outside the spin sector");
      if (!paths.full_ && !paths.covalent_ && paths.size_ == 0) return;
      beta_stack_.push_back({0, 0, 0, 0, true});
      ++counts_.prefixes;
      if (paths.words_ != 0) {
        for (int mode = 0; mode < beta_modes_; ++mode) {
          std::fill_n(beta_live_.data() + mode * paths.words_, paths.words_, ~Mask{0});
          if (paths.size_ % 64)
            beta_live_[mode * paths.words_ + paths.words_ - 1] =
                (Mask{1} << (paths.size_ % 64)) - 1;
        }
      }
      if (fixed) {
        for (int k = 0; k < paths.n_ && !beta_stack_.empty(); ++k)
          if (!extend_beta((beta >> k) & 1)) beta_stack_.clear();
      }
    }

    /** @brief Emit the next unique product and its depth-zero support label. */
    bool next(Mask& alpha, Mask& beta, bool& in_support) {
      for (;;) {
        if (alpha_active_) {
          if (next_alpha(alpha, beta, in_support)) return true;
          alpha_active_ = false;
          if (fixed_) return false;
          beta_stack_.pop_back(); // The completed beta leaf.
        }
        if (fixed_) {
          if (beta_stack_.empty() || fixed_started_) return false;
          fixed_started_ = true;
          start_alpha();
          continue;
        }
        if (!next_beta()) return false;
        start_alpha();
      }
    }

    Counts counts() const noexcept { return counts_; }

  private:
    struct Frame {
      int k, count, next;
      Mask mask;
      bool support;
    };

    bool extend_beta(int bit) {
      const Frame& parent = beta_stack_.back();
      const int k = parent.k, count = parent.count + bit;
      if (count > paths_->nb_ || count + paths_->n_ - k - 1 < paths_->nb_)
        return false;
      bool possible = true;
      if (paths_->covalent_) {
        // The earliest possible alpha string maximizes every alpha prefix.
        possible = count + std::min(k + 1, paths_->na_) + depth_ >= k + 1;
      } else if (!paths_->full_) {
        possible = false;
        const auto* src = beta_live_.data() + k * beta_modes_ * paths_->words_;
        auto* dst = beta_live_.data() + (k + 1) * beta_modes_ * paths_->words_;
        for (int mode = 0; mode < beta_modes_; ++mode) {
          const auto bound = paths_->offset(1, k, std::min(paths_->n_, count + mode));
          for (std::size_t w = 0; w < paths_->words_; ++w) {
            const Mask value = src[mode * paths_->words_ + w] &
                               paths_->bounds_[bound + w];
            dst[mode * paths_->words_ + w] = value;
            if (mode == beta_modes_ - 1) possible |= value != 0;
          }
          counts_.words += paths_->words_;
        }
      }
      if (!possible) return false;
      beta_stack_.push_back({k + 1, count, 0,
                             parent.mask | (static_cast<Mask>(bit) << k), true});
      ++counts_.prefixes;
      return true;
    }

    bool next_beta() {
      while (!beta_stack_.empty()) {
        Frame& top = beta_stack_.back();
        if (top.k == paths_->n_) return true;
        if (top.next == 2) {
          beta_stack_.pop_back();
          continue;
        }
        const int bit = top.next++;
        if (extend_beta(bit)) continue;
      }
      return false;
    }

    void start_alpha() {
      const int k = paths_->n_;
      if (paths_->words_ != 0) {
        const auto* beta = beta_live_.data() + k * beta_modes_ * paths_->words_;
        auto* alpha = alpha_live_.data();
        std::copy_n(beta, paths_->words_, alpha);
        if (depth_ == 1) {
          std::copy_n(beta, paths_->words_, alpha + paths_->words_);
          std::copy_n(beta + paths_->words_, paths_->words_,
                      alpha + 2 * paths_->words_);
        }
      }
      alpha_stack_.push_back({0, 0, 0, 0, true});
      alpha_active_ = true;
      ++counts_.prefixes;
    }

    bool extend_alpha(int bit) {
      const Frame& parent = alpha_stack_.back();
      const int k = parent.k, count = parent.count + bit;
      if (count > paths_->na_ || count + paths_->n_ - k - 1 < paths_->na_)
        return false;
      bool possible = true, support = true;
      if (paths_->covalent_) {
        const int beta_count = beta_stack_[k + 1].count;
        support = parent.support && count + beta_count >= k + 1;
        possible = count + beta_count + depth_ >= k + 1;
      } else if (!paths_->full_) {
        possible = false;
        support = false;
        const auto* src = alpha_live_.data() + k * alpha_modes_ * paths_->words_;
        auto* dst = alpha_live_.data() + (k + 1) * alpha_modes_ * paths_->words_;
        for (int mode = 0; mode < alpha_modes_; ++mode) {
          const auto bound = paths_->offset(0, k,
              std::min(paths_->n_, count + (mode == 1)));
          Mask present = 0;
          for (std::size_t w = 0; w < paths_->words_; ++w) {
            const Mask value = src[mode * paths_->words_ + w] &
                               paths_->bounds_[bound + w];
            dst[mode * paths_->words_ + w] = value;
            present |= value;
          }
          counts_.words += paths_->words_;
          possible |= present != 0;
          if (mode == 0) support = present != 0;
        }
      }
      if (!possible) return false;
      alpha_stack_.push_back({k + 1, count, 0,
                              parent.mask | (static_cast<Mask>(bit) << k), support});
      ++counts_.prefixes;
      return true;
    }

    bool next_alpha(Mask& alpha, Mask& beta, bool& in_support) {
      while (!alpha_stack_.empty()) {
        Frame& top = alpha_stack_.back();
        if (top.k == paths_->n_) {
          alpha = top.mask;
          beta = beta_stack_.back().mask;
          in_support = top.support;
          alpha_stack_.pop_back();
          ++counts_.products;
          counts_.support += in_support;
          return true;
        }
        if (top.next == 2) {
          alpha_stack_.pop_back();
          continue;
        }
        const int bit = top.next++;
        if (extend_alpha(bit)) continue;
      }
      return false;
    }

    const ProductPaths* paths_;
    int depth_;
    bool fixed_, fixed_started_ = false, alpha_active_ = false;
    int beta_modes_, alpha_modes_;
    Counts counts_;
    std::vector<Frame> beta_stack_, alpha_stack_;
    std::vector<Mask> beta_live_, alpha_live_;
  };

  Cursor cursor(int depth) const { return Cursor(*this, depth, false, 0); }
  Cursor column(Mask beta, int depth) const { return Cursor(*this, depth, true, beta); }

  template<class Emit>
  Counts visit(int depth, Emit&& emit) const {
    auto it = cursor(depth);
    Mask a, b;
    bool support;
    while (it.next(a, b, support)) emit(a, b, support);
    return it.counts();
  }

  template<class Emit>
  Counts visit_column(Mask beta, int depth, Emit&& emit) const {
    auto it = column(beta, depth);
    Mask a, b;
    bool support;
    while (it.next(a, b, support)) emit(a, b, support);
    return it.counts();
  }

  std::size_t metadata_bytes() const noexcept { return bounds_.capacity() * sizeof(Mask); }
  /** @brief Maximum bitset workspace, excluding the O(n) cursor frames. */
  std::size_t workspace_bytes(int depth) const noexcept {
    return (n_ + 1) * (depth == 0 ? 2 : 5) * words_ * sizeof(Mask);
  }

private:
  std::size_t offset(int spin, int k, int count) const noexcept {
    return ((static_cast<std::size_t>(spin) * n_ + k) * (n_ + 1) + count) * words_;
  }
  int n_, na_, nb_;
  bool covalent_ = false, full_ = false;
  std::size_t size_ = 0, words_ = 0;
  std::vector<Mask> bounds_;
};

} // namespace xmvb::mcscf
