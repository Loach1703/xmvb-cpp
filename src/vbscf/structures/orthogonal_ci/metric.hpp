#pragma once

#include <cstddef>
#include <vector>

#include <Eigen/Core>

namespace xmvb::vb {

class StructureAction;

/** @brief Exact complete-space overlap inverse with occupation-block factors.
 * Borrows an immutable action, which must outlive this object. Completeness of
 * the selected spin sector is required, not just completeness of its carrier.
 * No global dense structure overlap or expansion matrix is stored.
 */
class StructureMetricInverse {
public:
  StructureMetricInverse(
      const StructureAction& action,
      const std::vector<std::vector<int>>& alpha_strings,
      const std::vector<std::vector<int>>& beta_strings,
      bool complete_space,
      std::size_t factor_byte_limit,
      std::size_t workspace_byte_limit);
  StructureMetricInverse(const StructureMetricInverse&) = delete;
  StructureMetricInverse& operator=(const StructureMetricInverse&) = delete;
  StructureMetricInverse(StructureMetricInverse&&) noexcept = default;
  StructureMetricInverse& operator=(StructureMetricInverse&&) noexcept = default;

  /** @brief Applies the overlap inverse, chunking RHS by the workspace budget. */
  Eigen::MatrixXd apply(const Eigen::Ref<const Eigen::MatrixXd>& rhs) const;

  /** @brief Reuse the factors after their same-topology action has been moved.
   * No factor or integral data are copied. The new action must outlive this
   * inverse and retain the original spin-string/structure indexing.
   */
  void rebind(const StructureAction& action);
  std::size_t factor_bytes() const noexcept { return factor_bytes_; }
  std::size_t block_count() const noexcept { return blocks_.size(); }
  int max_block_size() const noexcept { return max_block_size_; }

  /** @brief Plans factor/workspace admission without allocating dense blocks. */
  static bool fits(
      const StructureAction& action,
      const std::vector<std::vector<int>>& alpha_strings,
      const std::vector<std::vector<int>>& beta_strings,
      std::size_t factor_byte_limit,
      std::size_t workspace_byte_limit);

private:
  struct Block {
    std::vector<int> indices;
    Eigen::MatrixXd lower;
  };

  static std::vector<std::vector<int>> occupation_groups(
      const StructureAction& action,
      const std::vector<std::vector<int>>& alpha_strings,
      const std::vector<std::vector<int>>& beta_strings);

  Eigen::MatrixXd solve_spin(const Eigen::Ref<const Eigen::MatrixXd>& rhs) const;
  const StructureAction* action_;
  std::vector<Block> blocks_;
  std::size_t factor_bytes_ = 0;
  int max_block_size_ = 0;
  int chunk_width_ = 0;
};

}  // namespace xmvb::vb
