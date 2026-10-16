#pragma once

#include "integrals/active/two_electron_result.hpp"
#include "integrals/ao/ri.hpp"
#include "vbscf/orbitals/preparation/result.hpp"

namespace xmvb::vb {

/**
 * @brief Transforms molecule-static AO-side RI factors into active-space RI factors.
 *
 * The RI factors are always retained for reverse-mode differentiation. When a
 * complete active-pair kernel and its GEMM workspace together fit within one
 * additional factor-matrix footprint, the result also carries that packed
 * lookup cache. Larger active spaces remain factorized and are consumed by
 * bounded active-pair tiles.
 */
class RiActiveSpaceTwoElectronBuilder {
public:
  xmvb::integrals::ActiveSpaceTwoElectronResult build(
      const xmvb::integrals::RiAoFactorization& ao_ri_result,
      const OrbitalPreparationResult& orbital_preparation_result,
      int n_basis_functions,
      int n_active_orbitals) const;
};

}  // namespace xmvb::vb
