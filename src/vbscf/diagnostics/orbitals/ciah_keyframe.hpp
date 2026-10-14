#pragma once

#include <functional>
#include <Eigen/Core>

namespace xmvb::vb {
class StructureAction;
class StructureTangentOperator;
struct AcceptedPointContext;
struct OrbitalPreparationInput;

namespace diagnostics {

/** @brief Opt-in first-keyframe derivative/path audit; never a solver policy. */
bool ciah_keyframe_audit_enabled();

/** @brief Frozen-integral carrier remainders at four fractions of the step. */
struct CiahCarrierAudit {
  Eigen::MatrixXd remainders;
};

/** @brief Measure carrier transport and its frozen-integral counterfactual. */
CiahCarrierAudit audit_ciah_orbital_geometry(
    const Eigen::VectorXd& orbital_displacement,
    const AcceptedPointContext& origin,
    const StructureAction& action,
    const StructureTangentOperator& structure,
    const std::function<OrbitalPreparationInput(const Eigen::VectorXd&)>& trial_orbitals);

/**
 * @brief Split finite structure response into orbital, CI and mixed terms.
 * Uses saved origin HVP images and matrix-free H/M actions, without dense
 * structure matrices. The callback evaluates a fixed-C gradient, not relaxed CI.
 */
void audit_ciah_keyframe(
    const Eigen::VectorXd& displacement,
    const Eigen::VectorXd& origin_gradient,
    const Eigen::VectorXd& hessian_displacement,
    const Eigen::VectorXd& actual_gradient,
    const Eigen::MatrixXd& origin_coefficients,
    double weight,
    const StructureTangentOperator& structure,
    const StructureAction& action,
    const CiahCarrierAudit& carrier,
    const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& evaluate_gradient);

}  // namespace diagnostics
}  // namespace xmvb::vb
