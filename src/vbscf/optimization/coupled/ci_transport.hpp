#pragma once

#include <memory>
#include <optional>
#include <Eigen/Core>

#include "vbscf/orbitals/charts/chart.hpp"
#include "vbscf/orbitals/preparation/overlap_pullback.hpp"
#include "vbscf/structures/overlap/representation.hpp"

namespace xmvb::vb {
struct AcceptedPointContext;
struct VbScfCoupledKeyframeEvaluation;
class StructureMetricInverse;
class StructureTangentOperator;
struct SelectedMetric;

/** @brief One selected-space retraction and its immutable trial geometry. */
struct CiTransportPoint {
  OrbitalPreparationInput input;
  OrbitalPreparationResult prepared;
  Eigen::MatrixXd upper;
  Eigen::MatrixXd origin_direction;
  Eigen::MatrixXd base_coefficients;
  Eigen::MatrixXd coefficients;
};

/** @brief Residual-weighted origin Hessian correction, split by input block. */
struct CiTransportCurvature {
  Eigen::VectorXd orbital_from_orbital;
  Eigen::VectorXd orbital_from_structure;
  Eigen::VectorXd structure_from_orbital;
};

/** @brief Projected finite CI transport; selected structures need not be closed. */
class CiTransport {
public:
  CiTransport(const OrbitalPreparationInput& input,
              const SparseParameterLayout& layout,
              const OrbitalChart& chart,
              const AcceptedPointContext& origin,
              const StructureAction& action,
              const StructureTangentOperator& structure);

  CiTransportPoint retract(const Eigen::VectorXd& orbital_step,
                           const Eigen::VectorXd& structure_step) const;
  Eigen::VectorXd gradient(const CiTransportPoint& point,
                          const VbScfCoupledKeyframeEvaluation& fixed) const;
  CiTransportCurvature curvature(const Eigen::VectorXd& orbital_direction,
                                const Eigen::VectorXd& structure_direction) const;
  Eigen::VectorXd origin_structure_gradient() const;

private:
  Eigen::MatrixXd solve(const Eigen::MatrixXd& upper,
                       const Eigen::MatrixXd& rhs,
                       const Eigen::MatrixXd* initial = nullptr) const;
  Eigen::VectorXd raw_direction(const Eigen::VectorXd& reduced) const;

  const OrbitalPreparationInput& input_;
  const SparseParameterLayout& layout_;
  const OrbitalChart& chart_;
  const AcceptedPointContext& origin_;
  const StructureAction& action_;
  const StructureTangentOperator& structure_;
  RepresentationOverlap overlap_;
  OrbitalOverlapPullback origin_pullback_;
  Eigen::MatrixXd upper_;
  Eigen::MatrixXd zero_;
  Eigen::MatrixXd c0_;
  Eigen::MatrixXd gc0_;
  Eigen::MatrixXd lambda0_;
  double coordinate_scale_;
  Eigen::VectorXd diagonal_;
  Eigen::MatrixXd residual_upper_gradient_;
  Eigen::MatrixXd residual_overlap_gradient_;
  const StructureMetricInverse* inverse_ = nullptr;
  /** One trial chart, shared by retraction and its gradient solve. */
  mutable Eigen::MatrixXd metric_upper_;
  mutable std::optional<RepresentationOverlap::PreparedMetric> prepared_metric_;
  mutable std::shared_ptr<const SelectedMetric> selected_metric_;
};

}  // namespace xmvb::vb
