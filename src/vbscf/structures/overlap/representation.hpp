#pragma once

#include <memory>
#include <vector>

#include <Eigen/Core>

namespace xmvb::vb {

class StructureAction;
class ExteriorOrbitalTransform;

/** @brief Chart cotangents actually consumed by the caller. */
enum class OverlapSide { Left, Right, Both };

/** Gradients in the independent upper-triangular entries of the two charts. */
struct UpperPair {
  Eigen::MatrixXd left;
  Eigen::MatrixXd right;
};

/** @brief Three origin-chart contractions sharing the same selected factors.
 * For f(L,R)=lambda^T F(L)^T F(R)c, directional is D_R grad_R f[dR],
 * mixed is grad_R f with c replaced by dc, and image is D F(R)^T[dR] F(R)lambda.
 */
struct OverlapTransportCurvature {
  Eigen::MatrixXd directional;
  Eigen::MatrixXd mixed;
  Eigen::VectorXd image;
};

/** @brief Finite-transport chart cotangents and reference CI image. */
struct OverlapTransportGradient {
  UpperPair charts;
  Eigen::VectorXd image;
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

  /** Adjoint of left_coeff^T F(left)^T F(right) right_coeff.
   * The unrequested side is empty, not computed and discarded.
   */
  UpperPair adjoint(
      const Eigen::Ref<const Eigen::MatrixXd>& left_upper,
      const Eigen::Ref<const Eigen::MatrixXd>& right_upper,
      const Eigen::Ref<const Eigen::VectorXd>& left_coeff,
      const Eigen::Ref<const Eigen::VectorXd>& right_coeff,
      OverlapSide side = OverlapSide::Both) const;

  /** Chart-direction derivative of adjoint, with coefficients held fixed. */
  UpperPair adjoint_direction(
      const Eigen::Ref<const Eigen::MatrixXd>& left_upper,
      const Eigen::Ref<const Eigen::MatrixXd>& right_upper,
      const Eigen::Ref<const Eigen::MatrixXd>& delta_left_upper,
      const Eigen::Ref<const Eigen::MatrixXd>& delta_right_upper,
      const Eigen::Ref<const Eigen::VectorXd>& left_coeff,
      const Eigen::Ref<const Eigen::VectorXd>& right_coeff,
      OverlapSide side = OverlapSide::Both) const;

  /** @brief Origin-chart directional/mixed cotangents and transpose action.
   * Uses a single bounded tile sweep for an orthogonal-CI carrier. The
   * selected structure columns need not be closed under orbital transforms.
   */
  OverlapTransportCurvature transport_curvature(
      const Eigen::Ref<const Eigen::MatrixXd>& upper,
      const Eigen::Ref<const Eigen::MatrixXd>& direction,
      const Eigen::Ref<const Eigen::VectorXd>& lambda,
      const Eigen::Ref<const Eigen::VectorXd>& coefficients,
      const Eigen::Ref<const Eigen::VectorXd>& coefficient_direction) const;

  /** @brief Pull back the finite selected-space transport in one tile sweep. */
  OverlapTransportGradient transport_gradient(
      const Eigen::Ref<const Eigen::MatrixXd>& current,
      const Eigen::Ref<const Eigen::MatrixXd>& origin,
      const Eigen::Ref<const Eigen::MatrixXd>& origin_direction,
      const Eigen::Ref<const Eigen::VectorXd>& lambda,
      const Eigen::Ref<const Eigen::VectorXd>& base,
      const Eigen::Ref<const Eigen::VectorXd>& reference,
      const Eigen::Ref<const Eigen::VectorXd>& transported) const;

private:
  const StructureAction* action_;
  const std::vector<std::vector<int>>& alpha_;
  const std::vector<std::vector<int>>& beta_;
  int n_active_;
  bool complete_carrier_;
};

}  // namespace xmvb::vb
