#pragma once

#include <cstddef>
#include <Eigen/Core>

namespace xmvb::integrals {

/** @brief Numerical workspace allowance for tiled RI orbital transforms. */
inline constexpr std::size_t ri_response_workspace_bytes=64ULL*1024*1024;

/** @brief Contract symmetric packed AO factors with one orbital index.
 * Output rows are (orbital, auxiliary), with auxiliary contiguous; columns
 * are AO functions. Workers own only one auxiliary-tile AO row slice.
 */
Eigen::MatrixXd transform_ri_one_index(
    const Eigen::Ref<const Eigen::MatrixXd>& packed_factors,
    const Eigen::Ref<const Eigen::MatrixXd>& coefficients,
    int first_aux,int count);

} // namespace xmvb::integrals
