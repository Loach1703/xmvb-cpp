#pragma once

#include <Eigen/Core>

namespace xmvb::dft {

/** Compute the restricted total density on one AO tile.
 * `projected` is reusable scratch with the same shape as `values`.
 */
void evaluate_density_tile(
    const Eigen::Ref<const Eigen::MatrixXd>& values,
    const Eigen::Ref<const Eigen::MatrixXd>& density_matrix,
    Eigen::MatrixXd* projected,
    Eigen::VectorXd* density);

}  // namespace xmvb::dft
