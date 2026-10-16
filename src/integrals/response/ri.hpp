#pragma once

#include <vector>

#include <Eigen/Core>

#include "integrals/active/two_electron_result.hpp"
#include "integrals/ao/ri.hpp"

namespace xmvb::integrals {

/**
 * @brief Accepted-point data for factor-native RI 2e responses.
 *
 * The cache borrows molecule-static whitened AO-pair factors, accepted
 * AO-by-active coefficients, and accepted active-pair factors `B = L T(C)`.
 * It owns only a bounded one-index transform, never AO-pair adjoints or an
 * orbital Hessian. Borrowed data must outlive every use of this cache.
 */
struct RiActiveTwoElectronResponseCache {
  int n_basis_functions = 0;
  int n_active_orbitals = 0;
  int n_auxiliary_functions = 0;
  std::vector<int> active_pair_first_indices;
  std::vector<int> active_pair_second_indices;
  const Eigen::MatrixXd* metric_whitened_ao_pair_factors = nullptr;
  const Eigen::MatrixXd* accepted_active_coefficients = nullptr;
  const Eigen::MatrixXd* accepted_active_pair_factors = nullptr;
  /** @brief Bounded one-index factors, row `a*n_aux+A`, column AO `mu`.
   * Empty means the identical transform is generated in auxiliary tiles.
   */
  Eigen::MatrixXd accepted_one_index_factors;
  /** @brief Resource-derived auxiliary tile width, unrelated to accuracy. */
  int auxiliary_tile_size = 0;
};

/**
 * @brief Builds a bounded one-index RI response cache with borrowed factors.
 *
 * @param ao_ri_factorization Molecule-static whitened AO-pair factors.
 * @param accepted_active_space_result Accepted RI active-space result carrying
 *        both `C` and `B`.
 * @param n_active_orbitals Number of active orbitals in `C`.
 */
RiActiveTwoElectronResponseCache build_ri_active_two_electron_response_cache(
    const RiAoFactorization& ao_ri_factorization,
    const ActiveSpaceTwoElectronResult& accepted_active_space_result,
    int n_active_orbitals);

/**
 * @brief Forms the active-pair factor direction `delta B = L delta T(C)[D]`.
 *
 * The returned matrix has auxiliary-function rows and packed active-pair
 * columns. No AO four-index integral tensor is formed or consulted.
 */
Eigen::MatrixXd compute_ri_active_pair_factor_directional_derivative(
    const RiActiveTwoElectronResponseCache& accepted_cache,
    const Eigen::Ref<const Eigen::MatrixXd>& dense_active_direction);

/**
 * @brief Forms several active-pair factor directions in one RI transform.
 *
 * All directions share the accepted `L_A C` transform. No AO-pair by
 * active-pair directional map is formed.
 */
std::vector<Eigen::MatrixXd>
compute_ri_active_pair_factor_directional_derivative_batch(
    const RiActiveTwoElectronResponseCache& accepted_cache,
    const std::vector<Eigen::MatrixXd>& dense_active_directions);

/**
 * @brief Forms packed `delta G = B^T delta B + delta B^T B`.
 *
 * @param accepted_cache Accepted-point cache containing `B`.
 * @param directional_active_pair_factors Factor direction `delta B` with the
 *        same shape as accepted `B`.
 */
std::vector<double>
compute_ri_packed_active_two_electron_integral_directional_derivative(
    const RiActiveTwoElectronResponseCache& accepted_cache,
    const Eigen::Ref<const Eigen::MatrixXd>&
        directional_active_pair_factors);

/**
 * @brief Forms packed `delta G` columns for a block of factor directions.
 */
Eigen::MatrixXd
compute_ri_packed_active_two_electron_integral_directional_derivative_batch(
    const RiActiveTwoElectronResponseCache& accepted_cache,
    const std::vector<Eigen::MatrixXd>& directional_active_pair_factors);

/**
 * @brief Pulls an arbitrary packed active-2e adjoint back to active orbitals.
 *
 * This general accepted-point pullback is suitable both for the ordinary
 * active-space adjoint and for an outer-response `delta` adjoint. It acts
 * directly through `B = L T(C)` and never reconstructs AO four-index ERIs.
 */
Eigen::MatrixXd backpropagate_ri_packed_active_two_electron_gradient(
    const std::vector<double>& packed_active_two_electron_gradient,
    const RiActiveTwoElectronResponseCache& accepted_cache);

/**
 * @brief Applies the fixed-packed-adjoint RI orbital Hessian to one direction.
 *
 * The result includes both the factor-response term driven by `delta B W`
 * and the derivative of the orbital pair-map pullback. A precomputed
 * `directional_active_pair_factors` may be supplied to share the factor action
 * with the packed-integral directional derivative.
 */
Eigen::MatrixXd
apply_ri_packed_active_two_electron_adjoint_hessian_vector(
    const std::vector<double>& packed_active_two_electron_gradient,
    const RiActiveTwoElectronResponseCache& accepted_cache,
    const Eigen::Ref<const Eigen::MatrixXd>& dense_active_direction,
    const Eigen::MatrixXd* directional_active_pair_factors = nullptr);

}  // namespace xmvb::integrals
