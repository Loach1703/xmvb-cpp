#pragma once

#include <vector>
#include <span>
#include <Eigen/Core>

#include "integrals/response/ri.hpp"

namespace xmvb::integrals {

/** @brief The two physical outputs required by the orbital Hessian pullback. */
struct RiAoResponse {
  Eigen::MatrixXd combined;  ///< G(2 D + (W C^T + C W^T)/2), symmetric covector.
  Eigen::MatrixXd projection; ///< G(D) C, without forming a second AO response.
};

/** @brief Shared auxiliary slabs and bounded narrow per-worker projections.
 * No worker owns an AO-square output or a copy of all RI factors.
 */
struct RiAoResponseWorkspace {
  struct Worker {
    Eigen::MatrixXd left, right, projection, core, mixed;
  };
  Eigen::MatrixXd direction, accepted;
  Eigen::VectorXd coulomb;
  std::vector<Worker> workers;
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

/** @brief Exact factor-native RI AO response for D=X A^T+A X^T.
 * @param active Borrowed accepted active coefficients and one-index factors.
 * @param inactive_duals A=C_i (C_i^T S C_i)^-1; may have zero columns.
 * @param horizontal_direction X=(I-P_i S) delta C_i.
 * @param weighted_active W=delta C_active (g_h+g_h^T).
 * @param workspace Reused shared slabs, charged to the RI response allowance.
 * No eigendecomposition, rank cutoff or strategy calibration is performed.
 */
RiAoResponse apply_ri_ao_response(
    const RiActiveTwoElectronResponseCache& active,
    const Eigen::Ref<const Eigen::MatrixXd>& inactive_duals,
    const Eigen::Ref<const Eigen::MatrixXd>& horizontal_direction,
    const Eigen::Ref<const Eigen::MatrixXd>& weighted_active,
    RiAoResponseWorkspace& workspace,
    std::size_t workspace_bytes=ri_response_workspace_bytes,
    RiInactiveFactorView inactive = {});

} // namespace xmvb::integrals
