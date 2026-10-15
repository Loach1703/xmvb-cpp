#include "vbscf/structures/assembly/action.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "mcscf/ci/factors.hpp"
#include "vbscf/structures/assembly/sources.hpp"

namespace xmvb::vb {

bool StructureAction::prefer_factor_response(const mcscf::DirectCiSigmaAction& sigma) const {
  const long double r = sources_->nodes.size();
  const long double l = r * (spin_exchange_parity_ ? 2 : 1);
  const long double u = static_cast<long double>(n_unique_alpha_) + n_unique_beta_;
  const long double d = sigma.channel_kernel().rows();
  const long double products = spin_exchange_parity_
      ? static_cast<long double>(n_unique_alpha_) * (n_unique_alpha_+spin_exchange_parity_)/2
      : static_cast<long double>(n_unique_alpha_) * n_unique_beta_;
  // Integral Gram: 2 U r d^2. Two residual generators: 8 U r d^2.
  // H/M source sweep: 2 U r d^2. Source-pair moments: 4 U l r (d+1)
  // for the integral adjoint and 8 U l r (d+2) for the two residuals.
  // Compare to one fused product sigma/Gram (4 products d^2).
  return r > 0 && 12*u*r*d*d + u*l*r*(12*d+20) < 4*products*d*d;
}

StructureResponseSweep StructureAction::source_response(
    const StructureAdjointState& state,
    const Eigen::Ref<const Eigen::MatrixXd>& coefficients,
    const ExteriorOrbitalTransform& alpha, const ExteriorOrbitalTransform& beta,
    const mcscf::DirectCiSigmaAction& sigma) const {
  const int n = sigma.n_orbitals(), d = sigma.channel_kernel().rows();
  const int n_states = state.weights.size();
  if (n_states <= 0 || state.structure_coefficients.rows() != n_structures_ ||
      state.structure_coefficients.cols() != n_states ||
      state.energies.size() != state.weights.size())
    throw std::invalid_argument("invalid selected-subspace adjoint state");
  StructureResponseSweep result{
      {Eigen::MatrixXd::Zero(n_structures_,n_states),
       Eigen::MatrixXd::Zero(n_structures_,n_states)},
      Eigen::MatrixXd::Zero(n,n),Eigen::MatrixXd::Zero(d,d),Eigen::MatrixXd::Zero(n,n)};
  for (int i = 0; i < n_states; ++i) {
    const double weight = state.weights[i], energy = state.energies[i];
    if (!std::isfinite(weight) || weight < 0 || !std::isfinite(energy))
      throw std::invalid_argument("invalid selected-subspace adjoint state");
    const Eigen::VectorXd response = expand_support(coefficients.col(i));
    if (response.isZero(0.0)) continue;
    // Expand original support once; images and all adjoints reuse these values.
    auto images = apply_source_columns(1,[&](int) { return response; },
        alpha,beta,sigma,false,mcscf::kFactorWorkspaceBytes);
    result.images.hamiltonian.col(i) = images.hamiltonian.col(0);
    result.images.overlap.col(i) = images.overlap.col(0);
    result.images.tiles.expansions += images.tiles.expansions;
    result.images.tiles.structure_projections += images.tiles.structure_projections;
    result.images.tiles.source_tiles += images.tiles.source_tiles;
    result.images.tiles.source_columns += images.tiles.source_columns;
    result.images.tiles.channel_blocks += images.tiles.channel_blocks;
    result.images.tiles.source_workspace_bound = std::max(
        result.images.tiles.source_workspace_bound,images.tiles.source_workspace_bound);
    if (weight == 0.0) continue;
    const Eigen::VectorXd accepted = expand_support(state.structure_coefficients.col(i));
    const auto factors = [&](const Eigen::VectorXd& values, bool physical) {
      return sources_->orthogonal_source(spin_products_,values,alpha,beta,physical);
    };
    const auto y = factors(accepted,false), physical_y = factors(accepted,true);
    const auto q = factors(response,false), physical_q = factors(response,true);
    const auto integral = mcscf::factor_integral_adjoint(sigma,physical_q,y);
    result.one_electron.noalias() += 2*weight*integral.one_electron;
    result.pair_kernel.noalias() += 2*weight*integral.pair_kernel;
    result.generator.noalias() += 2*weight*mcscf::factor_residual_generator_adjoint(
        sigma,physical_q,y,energy);
    result.generator.noalias() += 2*weight*mcscf::factor_residual_generator_adjoint(
        sigma,physical_y,q,energy);
  }
  return result;
}

} // namespace xmvb::vb
