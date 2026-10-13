#pragma once

#include <Eigen/Core>

namespace xmvb::vb {

/** @brief Normalized orbital tangent and its raw-coordinate left-inverse metric. */
struct WhitenedTangent {
  Eigen::MatrixXd basis;
  Eigen::MatrixXd inverse_gram;
};

/** @brief Whiten an Euclidean-orthonormal quotient, retaining its horizontal gauge.
 * The normalization derivative supplies a rank-one correction. Its raw Gram
 * inverse is recovered from the same spectrum and a two-dimensional update,
 * without assembling or factoring another reduced-size matrix.
 */
WhitenedTangent whiten_orbital_tangent(
    const Eigen::Ref<const Eigen::MatrixXd>& quotient,
    const Eigen::Ref<const Eigen::VectorXd>& coefficients,
    const Eigen::Ref<const Eigen::MatrixXd>& overlap,
    bool sphere,
    const Eigen::MatrixXd* quotient_metric = nullptr);

}  // namespace xmvb::vb
