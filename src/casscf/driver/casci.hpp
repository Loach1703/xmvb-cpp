#pragma once

#include <cstddef>
#include <vector>
#include "casscf/integrals/active.hpp"
#include "linalg/eigensolver.hpp"

namespace xmvb::casscf {

/** @brief Fixed-orbital CI accuracy and total retained solver-workspace budget. */
struct CasciOptions {
  int n_roots = 1;
  double energy_tolerance = 1e-10;
  double residual_tolerance = 1e-10;
  double stationarity_tolerance = 1e-8;
  std::size_t workspace_bytes = 256ULL*1024*1024;
};

struct CasciResult {
  /** Active-space eigenpairs; images and eigenvalues exclude core_energy. */
  linalg::DavidsonResult active_ci;
  double core_energy = 0.0;
  std::vector<double> spin_squared;
  /** @brief Total Born--Oppenheimer energy, including inactive and nuclear terms. */
  double energy(int root = 0) const {
    return core_energy + active_ci.eigenpairs.eigenvalues.at(root);
  }
};

/** @brief Pure-spin CASCI using the common direct-CI action and Davidson.
 * twice_spin is 2S. The carrier uses Ms=S; guesses and every preconditioned
 * direction are projected to S, without a spin penalty or any VB structures.
 * Coefficients are column-major alpha/beta products, packed one root per column.
 */
CasciResult solve_casci(const CasIntegrals& integrals, int n_electrons,
    int twice_spin, const CasciOptions& options = {});

}  // namespace xmvb::casscf
