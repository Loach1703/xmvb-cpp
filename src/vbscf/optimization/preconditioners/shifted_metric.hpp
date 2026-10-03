#pragma once

#include <functional>

#include <Eigen/Core>

namespace xmvb::vb {

using PreconditionerAction =
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;

/**
 * @brief Approximately inverts the shifted local orbital model for NEO.
 *
 * The caller supplies a dimensionless residual fraction relative to the
 * current orbital covector. The inner solve measures both sides in the dual
 * norm of the positive block inverse. The outer NEO solver checks the complete
 * physical coupled KKT residual and remains responsible for step accuracy.
 */
Eigen::VectorXd apply_shifted_metric_preconditioner(
    const Eigen::VectorXd& covector,
    double shift,
    double relative_residual_target,
    const PreconditionerAction& apply_model,
    const PreconditionerAction& apply_metric,
    const PreconditionerAction& apply_block_inverse);

}  // namespace xmvb::vb
