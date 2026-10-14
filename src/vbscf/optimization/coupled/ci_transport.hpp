#pragma once

#include <memory>
#include <Eigen/Core>

#include "vbscf/orbitals/charts/chart.hpp"
#include "vbscf/orbitals/preparation/overlap_pullback.hpp"
#include "vbscf/structures/overlap/representation.hpp"

namespace xmvb::vb {
struct AcceptedPointContext;
struct VbScfCoupledKeyframeEvaluation;
class StructureMetricInverse;
class StructureTangentOperator;

/** @brief One selected-space retraction and its immutable trial geometry. */
struct CiTransportPoint {
  OrbitalPreparationInput input;
  OrbitalPreparationResult prepared;
  Eigen::MatrixXd upper;
  Eigen::MatrixXd origin_direction;
  Eigen::VectorXd base_coefficients;
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
  Eigen::VectorXd solve(const Eigen::MatrixXd& upper,
                       const Eigen::VectorXd& rhs) const;
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
  Eigen::VectorXd c0_;
  Eigen::VectorXd gc0_;
  Eigen::VectorXd lambda0_;
  Eigen::VectorXd diagonal_;
  Eigen::MatrixXd residual_upper_gradient_;
  Eigen::MatrixXd residual_overlap_gradient_;
  const StructureMetricInverse* inverse_ = nullptr;
};

}  // namespace xmvb::vb
