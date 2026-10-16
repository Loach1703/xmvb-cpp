#include "vbscf/structures/evaluation/evaluator.hpp"

namespace xmvb::vb {

StructureAccumulationResult StructureMatrixEvaluator::evaluate(
    const VbScfInput& input,
    const PreparedActiveSpaceContext& prepared_active_space) const {
  return structure_builder_.build(
      input.structure_data.alpha_det,
      input.structure_data.beta_det,
      input.structure_data.determinant_to_structure_terms,
      prepared_active_space.orbital_result.active_orbital_overlap_matrix,
      prepared_active_space.active_space_one_electron_result.h1e_act,
      input.orbital_preparation_input.n_active_orbitals,
      prepared_active_space.active_space_two_electron_result,
      input.structure_data.n_structures);
}

SameSpinPairCacheContext StructureMatrixEvaluator::build_pair_cache(
    const VbScfInput& input,
    const PreparedActiveSpaceContext& prepared_active_space,
    SameSpinPairCacheBuildOptions options) const {
  return build_same_spin_pair_cache_context(
      input.structure_data.alpha_det,
      input.structure_data.beta_det,
      structure_builder_.make_pair_evaluator(),
      prepared_active_space.orbital_result.active_orbital_overlap_matrix,
      prepared_active_space.active_space_one_electron_result.h1e_act,
      input.orbital_preparation_input.n_active_orbitals,
      prepared_active_space.active_space_two_electron_result,
      options);
}

}  // namespace xmvb::vb
