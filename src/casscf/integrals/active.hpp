#pragma once

#include <functional>
#include "mcscf/ci/integrals.hpp"

namespace xmvb::casscf {

/** @brief Active Hamiltonian and the inactive-electron plus nuclear constant. */
struct CasIntegrals {
  mcscf::ActiveIntegrals active;
  double core_energy = 0.0;
};

/** @brief Read one chemists' MO Coulomb integral without requiring an MO^4 array. */
using MoCoulomb = std::function<double(int,int,int,int)>;

/** @brief Restrict orthonormal MO integrals to a complete active space.
 * one_electron and coulomb refer to the full MO basis, with doubly
 * occupied inactive orbitals first and active orbitals immediately after.
 * The returned active one-electron operator includes the inactive potential.
 * Only requested core/active entries are read; virtual four-index integrals
 * are neither requested nor stored.
 */
CasIntegrals make_active_integrals(
    const Eigen::Ref<const Eigen::MatrixXd>& one_electron,
    const MoCoulomb& coulomb, int n_inactive, int n_active,
    double nuclear_energy);

}  // namespace xmvb::casscf
