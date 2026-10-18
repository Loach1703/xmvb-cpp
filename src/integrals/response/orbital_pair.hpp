#pragma once

#include <span>
#include <Eigen/Core>

namespace xmvb::integrals {

using OrbitalPairMatrix =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

/** @brief Pull back an integral adjoint through shared (AO active|active pair).
 * Products have rows (AO, active) and packed-pair columns. The symmetric pair
 * adjoint doubles its diagonal; the orbital-pair derivative supplies the two
 * remaining index placements, including coincident indices.
 */
Eigen::MatrixXd pullback_orbital_pairs(
    const OrbitalPairMatrix& products,
    const OrbitalPairMatrix& pair_adjoint,
    std::span<const int> first,
    std::span<const int> second,
    int n_active);

}  // namespace xmvb::integrals
