#pragma once

#include <memory>
#include <vector>

#include <Eigen/Core>

namespace xmvb::vb {

class StructureAction;
class ExteriorOrbitalTransform;

/** Gradients in the independent upper-triangular entries of the two charts. */
struct UpperPair {
  Eigen::MatrixXd left;
  Eigen::MatrixXd right;
};

/** Cross Gram action F(R_left)^T F(R_right), with F(R) = wedge(R) W.
 *
 * W is the immutable selected-structure expansion retained by `action`. The
 * selected structure space may be arbitrary. The action and both supplied
 * unique-string lists must outlive this object; their ordering must remain
 * unchanged throughout that lifetime.
 */
class RepresentationOverlap {
public:
  RepresentationOverlap(
      const StructureAction& action,
      const std::vector<std::vector<int>>& alpha_strings,
      const std::vector<std::vector<int>>& beta_strings,
      int n_active);

  /** Reusable symmetric Gram action for one fixed orbital chart.
   * The action and spin-string lists supplied to RepresentationOverlap must
   * outlive this object. A selected Gram fitting the existing workspace is
   * projected once; larger selections retain the exact factor-tile action.
   */
  class PreparedMetric {
  public:
    Eigen::MatrixXd operator()(
        const Eigen::Ref<const Eigen::MatrixXd>& structure_coeff_block) const;
    /** @brief Bounded selected Gram, empty when the operator is streamed. */
    const Eigen::MatrixXd& selected_matrix() const noexcept { return selected_gram_; }
    /** @brief Whether the selected matrix fits its existing workspace quarter. */
    static bool fits_selected(int dimension) noexcept;

  private:
    friend class RepresentationOverlap;
    PreparedMetric(const StructureAction& action,
                   const std::vector<std::vector<int>>& alpha_strings,
                   const std::vector<std::vector<int>>& beta_strings,
                   bool complete_carrier,
                   const Eigen::Ref<const Eigen::MatrixXd>& upper);

    const StructureAction* action_;
    const std::vector<std::vector<int>>* alpha_;
    const std::vector<std::vector<int>>* beta_;
    std::shared_ptr<const ExteriorOrbitalTransform> alpha_transform_;
    std::shared_ptr<const ExteriorOrbitalTransform> beta_transform_;
    Eigen::MatrixXd kernel_;
    Eigen::MatrixXd selected_gram_;
  };

  PreparedMetric prepare_metric(
      const Eigen::Ref<const Eigen::MatrixXd>& upper) const;

  Eigen::MatrixXd apply(
      const Eigen::Ref<const Eigen::MatrixXd>& left_upper,
      const Eigen::Ref<const Eigen::MatrixXd>& right_upper,
      const Eigen::Ref<const Eigen::MatrixXd>& structure_coeff_block) const;

  /** Chart-direction derivative of the cross Gram action. */
  Eigen::MatrixXd apply_direction(
      const Eigen::Ref<const Eigen::MatrixXd>& left_upper,
      const Eigen::Ref<const Eigen::MatrixXd>& right_upper,
      const Eigen::Ref<const Eigen::MatrixXd>& delta_left_upper,
      const Eigen::Ref<const Eigen::MatrixXd>& delta_right_upper,
      const Eigen::Ref<const Eigen::MatrixXd>& structure_coeff_block) const;

  /** Adjoint of left_coeff^T F(left)^T F(right) right_coeff. */
  UpperPair adjoint(
      const Eigen::Ref<const Eigen::MatrixXd>& left_upper,
      const Eigen::Ref<const Eigen::MatrixXd>& right_upper,
      const Eigen::Ref<const Eigen::VectorXd>& left_coeff,
      const Eigen::Ref<const Eigen::VectorXd>& right_coeff) const;

  /** Chart-direction derivative of adjoint, with coefficients held fixed. */
  UpperPair adjoint_direction(
      const Eigen::Ref<const Eigen::MatrixXd>& left_upper,
      const Eigen::Ref<const Eigen::MatrixXd>& right_upper,
      const Eigen::Ref<const Eigen::MatrixXd>& delta_left_upper,
      const Eigen::Ref<const Eigen::MatrixXd>& delta_right_upper,
      const Eigen::Ref<const Eigen::VectorXd>& left_coeff,
      const Eigen::Ref<const Eigen::VectorXd>& right_coeff) const;

private:
  const StructureAction* action_;
  const std::vector<std::vector<int>>& alpha_;
  const std::vector<std::vector<int>>& beta_;
  int n_active_;
  bool complete_carrier_;
};

}  // namespace xmvb::vb
