#pragma once

#include <functional>
#include <Eigen/Core>

namespace xmvb::vb {
struct AcceptedPointContext;
struct OrbitalPreparationInput;
struct FullDeterminantStructureData;
struct VbScfCoupledKeyframeEvaluation;
class SparseParameterLayout;
class StructureAction;
class StructureTangentOperator;

namespace diagnostics {

/** @brief Select the finite-transport diagnostic, never a solver policy. */
bool ciah_ci_transport_audit_enabled();

/** @brief Compare a transported first keyframe without changing its candidate. */
void audit_ciah_ci_transport(
    const Eigen::VectorXd& displacement,
    const Eigen::VectorXd& predicted_structure,
    const Eigen::VectorXd& actual_structure,
    const AcceptedPointContext& origin,
    const FullDeterminantStructureData& structures,
    const StructureAction& action,
    const StructureTangentOperator& structure,
    const SparseParameterLayout& layout,
    const std::function<OrbitalPreparationInput(const Eigen::VectorXd&)>& trial_orbitals,
    const std::function<VbScfCoupledKeyframeEvaluation(
        const Eigen::VectorXd&, const Eigen::MatrixXd&)>& evaluate);

}  // namespace diagnostics
}  // namespace xmvb::vb
