#pragma once

#include <cstddef>
#include <Eigen/Core>
#include "core/workspace.hpp"

namespace xmvb::integrals {

/** @brief Contract symmetric packed AO factors with one orbital index.
 * Output rows are (orbital, auxiliary), with auxiliary contiguous; columns
 * are AO functions. Workers own only one auxiliary-tile AO row slice.
 */
Eigen::MatrixXd transform_ri_one_index(
    const Eigen::Ref<const Eigen::MatrixXd>& packed_factors,
    const Eigen::Ref<const Eigen::MatrixXd>& coefficients,
    int first_aux,int count);

/** @brief Fill a resident image with bounded transform scratch.
 * The caller admits the destination to its accepted-point storage budget.
 */
void fill_ri_one_index_tiled(
    const Eigen::Ref<const Eigen::MatrixXd>& packed_factors,
    const Eigen::Ref<const Eigen::MatrixXd>& coefficients,
    Eigen::Ref<Eigen::MatrixXd> result,
    std::size_t workspace_bytes = workspace::available());

} // namespace xmvb::integrals
