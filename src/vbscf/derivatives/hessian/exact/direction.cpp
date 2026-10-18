#include "vbscf/derivatives/hessian/exact/apply_internal.hpp"
#include "vbscf/derivatives/hessian/exact/state_internal.hpp"

#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "vbscf/derivatives/hessian/responses/active_space/integral_direction.hpp"
#include "vbscf/derivatives/hessian/responses/orbital/preparation.hpp"
#include "integrals/response/ri.hpp"

namespace xmvb::vb {

ExactHvpOperator::State::PrecomputedDirection
ExactHvpOperator::State::prepare_direction(
    Eigen::VectorXd packed_direction) const {
  const auto& orbital_input = current_input_->orbital_preparation_input;
  const int n_inactive =
      (orbital_input.n_total_electrons - orbital_input.n_active_electrons) / 2;
  PrecomputedDirection result;
  result.packed_direction = std::move(packed_direction);
  result.dense_orbital_tangent_context =
      build_dense_orbital_tangent_context(
          orbital_input,
          parameter_view_,
          result.packed_direction,
          *accepted_orbital_preparation_cache_);
  result.orbital_preparation_directional_result =
      build_orbital_preparation_directional_result(
          orbital_input,
          result.dense_orbital_tangent_context,
          n_inactive,
          orbital_input.n_active_orbitals,
          *accepted_orbital_preparation_cache_,accepted_ri_factorization_==nullptr);
  return result;
}

integrals::RiAoResponse ExactHvpOperator::State::ri_ao_response(
    const PrecomputedDirection& direction,
    const Eigen::Ref<const Eigen::MatrixXd>& weighted_active) const {
  const auto& geometry=*accepted_orbital_preparation_cache_;
  const auto& input=current_input_->orbital_preparation_input;
  const int n=input.n_basis_functions;
  const int ni=(input.n_total_electrons-input.n_active_electrons)/2;
  // Without an inactive sector, neither its density derivative nor its
  // cotangent is consumed by the physical orbital pullback.
  if (!ni) return {Eigen::MatrixXd::Zero(n,n),
      Eigen::MatrixXd::Zero(n,input.n_active_orbitals)};
  const auto dci=direction.dense_orbital_tangent_context.delta_normalized_orbitals.leftCols(ni);
  const auto sci=geometry.basis_overlap_times_normalized.leftCols(ni);
  std::vector<int> moving;
  for (int i=0;i<ni;++i) if (!dci.col(i).isZero(0)) moving.push_back(i);
  Eigen::MatrixXd delta(n,moving.size()),duals(n,moving.size());
  for (int i=0;i<static_cast<int>(moving.size());++i) {
    delta.col(i)=dci.col(moving[i]);
    duals.col(i)=geometry.inactive_auxiliary.col(moving[i]);
  }
  // A column absent from delta C_i contributes exactly zero to delta P_i.
  // This is structural zero pruning, with no numerical rank threshold.
  const Eigen::MatrixXd horizontal=delta-geometry.inactive_auxiliary*(sci.transpose()*delta);
  return integrals::apply_ri_ao_response(*accepted_ri_two_electron_cache_,
      duals,horizontal,weighted_active,ri_ao_response_workspace_,
      integrals::ri_response_workspace_bytes,
      {accepted_ri_inactive_one_index_factors_.size() != 0
          ? &accepted_ri_inactive_one_index_factors_ : nullptr, moving});
}

ActiveSpaceIntegralDirectionView
ExactHvpOperator::State::build_active_integral_direction(
    const OrbitalPreparationDirectionalResult& orbital_direction,
    const Eigen::Ref<const Eigen::MatrixXd>& delta_h1e_times_active,
    const Eigen::VectorXd* precomputed_delta_packed_two_electron,
    ActiveSpaceIntegralDirectionWorkspace* workspace,
    bool packed_two_electron_required) const {
  const auto& accepted_two_electron = accepted_point_context_
      ->prepared_active_space.active_space_two_electron_result;
  const ActiveSpaceIntegralDirectionContext integral_context{
      *current_input_,
      accepted_ri_two_electron_cache_.has_value()
          ? nullptr
          : &accepted_exact_two_electron_cache_,
      accepted_active_auxiliary_orbitals_,
      accepted_basis_overlap_times_active_auxiliary_orbitals_,
      accepted_ao_effective_one_electron_times_active_auxiliary_orbitals_,
      accepted_ao_effective_one_electron_transpose_times_active_auxiliary_orbitals_,
      accepted_two_electron.dense_active_coefficients};
  return xmvb::vb::build_active_space_integral_direction(
      integral_context,
      ActiveSpaceIntegralTangent{
          orbital_direction.delta_active_auxiliary_orbitals,
          orbital_direction.delta_active_auxiliary_orbitals,
          delta_h1e_times_active,
          precomputed_delta_packed_two_electron,
          packed_two_electron_required},
      workspace);
}

}  // namespace xmvb::vb
