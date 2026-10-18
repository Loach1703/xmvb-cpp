#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include <Eigen/Core>

namespace xmvb::vb {

/** @brief Exact compound-matrix rows by descending occupied-row recursion.
 * Only Pascal coefficients and the caller's column permutation are retained.
 * Each worker uses two single-spin vectors, never a spin-product matrix.
 */
class ExteriorRowPlan {
public:
  ExteriorRowPlan(int n_orbitals, int n_electrons,
      std::span<const std::uint64_t> masks);

  Eigen::MatrixXd apply(const Eigen::MatrixXd& upper,
      std::span<const std::uint64_t> masks, std::span<const int> rows) const;
  std::size_t dynamic_bytes() const noexcept;

private:
  std::uint64_t choose(int n, int k) const noexcept;
  int n_, electrons_;
  std::vector<std::uint64_t> binomial_;
  std::vector<int> column_ranks_;
};

} // namespace xmvb::vb
