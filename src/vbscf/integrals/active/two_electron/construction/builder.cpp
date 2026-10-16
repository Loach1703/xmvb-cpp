#include "vbscf/integrals/active/two_electron/construction/builder.hpp"

#include "integrals/active/exact.hpp"

namespace xmvb::vb {

xmvb::integrals::ActiveSpaceTwoElectronResult ActiveSpaceTwoElectronBuilder::build(
    const xmvb::integrals::AoIntegralInput& ao,
    const OrbitalPreparationResult& orbitals, int n_active) const {
  return integrals::transform_exact_active_integrals(
      ao, orbitals.active_sparse_row_offsets,
      orbitals.active_sparse_orbital_indices,
      orbitals.active_sparse_values, n_active);
}

}  // namespace xmvb::vb
