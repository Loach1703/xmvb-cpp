#pragma once

#include <cstddef>
#include <vector>

#include <Eigen/Core>

namespace xmvb::vb {

class StructureAction;

/** @brief Orbital-independent spin Gram G=W^T W, factored by occupation.
 * Valid for any independent selection of VB structures with a definite spatial
 * occupation per structure; neither CAS completeness nor orbital closure is
 * assumed. Stores sum of squared occupation-block sizes, not a global Gram.
 */
class SpinMetric {
public:
  enum class Part { Gram, Lower, Upper };

  SpinMetric(const StructureAction& action,
      const std::vector<std::vector<int>>& alpha,
      const std::vector<std::vector<int>>& beta, std::size_t byte_limit);

  /** @brief Exact G, L or L^T solve for G=L L^T. No orbital transform. */
  Eigen::MatrixXd solve(const Eigen::Ref<const Eigen::MatrixXd>& rhs,
      Part part = Part::Gram) const;

  /** @brief SPD block-scaled inverse, matching the supplied diagonal sums.
   * One positive scalar per occupation preserves the full spin Gram geometry.
   * The supplied diagonal may be the existing positive Jacobi estimate of M.
   * This is only a preconditioner of the physical nonorthogonal metric M.
   * Its use never replaces M, changes the selected space or the solve target.
   */
  Eigen::VectorXd precondition(const Eigen::VectorXd& rhs,
      const Eigen::VectorXd& physical_diagonal) const;

  std::size_t factor_bytes() const noexcept { return factor_bytes_; }
  std::size_t block_count() const noexcept { return blocks_.size(); }
  int max_block_size() const noexcept { return max_block_size_; }
  const Eigen::VectorXd& diagonal() const noexcept { return diagonal_; }

  /** @brief Factor residency only, never admission to the physical CI path. */
  static bool fits(const StructureAction& action,
      const std::vector<std::vector<int>>& alpha,
      const std::vector<std::vector<int>>& beta, std::size_t byte_limit);

private:
  struct Block {
    std::vector<int> indices;
    Eigen::MatrixXd lower;
  };
  static std::vector<std::vector<int>> occupation_groups(
      const StructureAction& action,
      const std::vector<std::vector<int>>& alpha,
      const std::vector<std::vector<int>>& beta);

  int dimension_;
  Eigen::VectorXd diagonal_;
  std::vector<Block> blocks_;
  std::size_t factor_bytes_ = 0;
  int max_block_size_ = 0;
};

}  // namespace xmvb::vb
