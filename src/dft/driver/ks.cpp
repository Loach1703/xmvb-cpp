#include "dft/driver/ks.hpp"

#include <span>
#include <stdexcept>

#include "dft/integration/potential.hpp"
#include "dft/xc/functional.hpp"
#include "integrals/jk/contraction.hpp"

namespace xmvb::dft {

KsResult solve_ks(int n_alpha, int n_beta,
    const integrals::LibcintInput& basis,
    const Eigen::Ref<const Eigen::MatrixXd>& overlap,
    const integrals::AoIntegralInput& ao,
    const KsOptions& options) {
  if (!options.unrestricted && n_alpha != n_beta)
    throw std::invalid_argument("RKS requires equal alpha and beta populations");
  const int n_spin = options.unrestricted ? 2 : 1;
  const Functional functional(options.functional, n_spin);
  const auto grid = make_grid(basis, options.grid);
  const XcIntegrator integrator(basis, grid, functional);
  const auto& core = ao.ao_core_hamiltonian_matrix;
  KsResult result;
  result.grid_points = static_cast<std::size_t>(grid.weights.size());
  bool first = true;
  const auto retain = [&](double energy, const XcResult& xc) {
    if (first) {
      result.initial_electronic_energy = energy;
      first = false;
    }
    result.xc_energy = xc.energy;
    result.integrated_electrons = xc.electron_count;
  };
  const auto energy = [&](const Eigen::Ref<const Eigen::MatrixXd>& total,
                          const Eigen::MatrixXd& coulomb, double exc) {
    return (total.array() * core.array()).sum() +
        0.5 * (total.array() * coulomb.array()).sum() + exc;
  };

  if (options.unrestricted) {
    const auto evaluate = [&](const scf::SpinMatrices& density) {
      const std::array<Eigen::Ref<const Eigen::MatrixXd>, 2> channels{
          density[0], density[1]};
      const auto xc = integrator.evaluate(channels);
      const Eigen::MatrixXd total = density[0] + density[1];
      const Eigen::MatrixXd coulomb = integrals::contract_jk(total, ao);
      scf::UnrestrictedEvaluation physical;
      physical.electronic_energy = energy(total, coulomb, xc.energy);
      for (int spin = 0; spin < 2; ++spin)
        physical.fock[spin] = core + coulomb + xc.potential[spin];
      retain(physical.electronic_energy, xc);
      return physical;
    };
    const scf::ApplyUnrestrictedResponse response=[&](const scf::SpinMatrices& density,
        const scf::SpinMatrices& direction) {
      const std::array<Eigen::Ref<const Eigen::MatrixXd>,2> reference{density[0],density[1]};
      const std::array<Eigen::Ref<const Eigen::MatrixXd>,2> delta{direction[0],direction[1]};
      auto result=integrator.response(reference,delta);
      const Eigen::MatrixXd total=direction[0]+direction[1];
      const auto j=integrals::contract_jk(total,ao);
      for (auto& spin:result) spin+=j;
      return result;
    };
    result.state = scf::solve_unrestricted(n_alpha, n_beta, overlap, core,
        evaluate, options.scf,{},response);
  } else {
    const auto evaluate = [&](const Eigen::Ref<const Eigen::MatrixXd>& density) {
      const std::array<Eigen::Ref<const Eigen::MatrixXd>, 1> channels{density};
      const auto xc = integrator.evaluate(channels);
      const Eigen::MatrixXd coulomb = integrals::contract_jk(density, ao);
      scf::Evaluation physical;
      physical.electronic_energy = energy(density, coulomb, xc.energy);
      physical.fock = core + coulomb + xc.potential[0];
      retain(physical.electronic_energy, xc);
      return physical;
    };
    const scf::ApplyResponse response=[&](const Eigen::Ref<const Eigen::MatrixXd>& density,
        const Eigen::Ref<const Eigen::MatrixXd>& direction) {
      const std::array<Eigen::Ref<const Eigen::MatrixXd>,1> reference{density};
      const std::array<Eigen::Ref<const Eigen::MatrixXd>,1> delta{direction};
      auto result=integrator.response(reference,delta);
      result[0]+=integrals::contract_jk(direction,ao);
      return std::move(result[0]);
    };
    result.state = scf::solve_restricted(n_alpha + n_beta, overlap, core,
        evaluate, options.scf,{},response);
  }
  return result;
}

}  // namespace xmvb::dft
