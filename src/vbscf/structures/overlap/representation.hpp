#pragma once

#include <vector>

#include <Eigen/Core>

namespace xmvb::vb {

class StructureAction;

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
