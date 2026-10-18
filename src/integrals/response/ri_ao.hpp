#pragma once

#include <span>
#include <vector>
#include <Eigen/Core>

#include "integrals/response/ri.hpp"

namespace xmvb::integrals {

/** @brief Narrow AO responses needed by the orbital Hessian pullback. */
struct RiAoProjectedResponse {
  Eigen::MatrixXd inactive; ///< 2 G(2 D + K/2) A0, AO by all inactive.
  Eigen::MatrixXd active;   ///< G(D) C, AO by active.
  Eigen::MatrixXd factor_direction{}; ///< Optional delta B, auxiliary by active pair.
  Eigen::MatrixXd fixed_two_electron{}; ///< Optional fixed-adjoint image, AO by active.
};

/** @brief Borrowed fixed-2e request in the accepted active coefficient frame.
 * `weight` maps the active direction to the AO-response density weight:
 * Wbar = direction * weight. `adjoint` is the accepted packed 2e cotangent.
 */
struct RiActiveHvpInput {
  const Eigen::MatrixXd& direction;
  const Eigen::MatrixXd& weight;
  const std::vector<double>& adjoint;
};

/** @brief Shared auxiliary slabs and narrow tile projections.
 * AO projections have disjoint writers. Only the fixed-jet reduction uses
 * small AO-by-active worker accumulators; integral slabs remain shared.
 */
struct RiAoResponseWorkspace {
  Eigen::MatrixXd direction, accepted;
  Eigen::MatrixXd direction_projection, accepted_projection, charges;
  Eigen::VectorXd coulomb;
};

/** @brief Borrowed accepted inactive transform, with exact column selection. */
struct RiInactiveFactorView {
  const Eigen::MatrixXd* factors = nullptr;
  std::span<const int> orbitals;
};

/** @brief G(C W C^T) from accepted L_A C and active-pair factors.
 * W must be symmetric. No AO eigendecomposition or repeated L_A C transform
 * is needed when the accepted factors are resident; otherwise the same
 * operator consumes bounded auxiliary tiles.
 */
Eigen::MatrixXd apply_ri_active_density_response(
    const RiActiveTwoElectronResponseCache& active,
    const Eigen::Ref<const Eigen::MatrixXd>& weights);

/** @brief Compute only the AO-by-orbital projections consumed by the HVP.
 * `moving_inactive_duals` and `horizontal_direction` contain matching moving
 * columns; `full_inactive_duals` contains every accepted inactive dual A0.
 * When `inactive.factors` is present it stores L_A A0 as rows
 * `orbital*n_auxiliary+auxiliary`, and `inactive.orbitals` maps each moving
 * column to its A0 column. Without cached factors, only the moving duals
 * are transformed; a full L_A A0 transform is never required. The returned
 * inactive projection includes the factor two required by the symmetric
 * inactive-density pullback. No AO-square response is materialized.
 * With `fixed_hvp`, `weighted_active` must equal
 * `fixed_hvp->direction * fixed_hvp->weight`; the returned optional images
 * are in the accepted active-coefficient frame.
 */
RiAoProjectedResponse apply_ri_projected_ao_response(
    const RiActiveTwoElectronResponseCache& active,
    const Eigen::Ref<const Eigen::MatrixXd>& moving_inactive_duals,
    const Eigen::Ref<const Eigen::MatrixXd>& horizontal_direction,
    const Eigen::Ref<const Eigen::MatrixXd>& weighted_active,
    const Eigen::Ref<const Eigen::MatrixXd>& full_inactive_duals,
    RiAoResponseWorkspace& workspace,
    std::size_t workspace_bytes = ri_response_workspace_bytes,
    RiInactiveFactorView inactive = {},
    const RiActiveHvpInput* fixed_hvp = nullptr);

} // namespace xmvb::integrals
