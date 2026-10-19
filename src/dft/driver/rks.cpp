#include "dft/driver/rks.hpp"

#include "dft/integration/potential.hpp"
#include "dft/xc/functional.hpp"
#include "integrals/jk/contraction.hpp"

namespace xmvb::dft {

RksResult solve_rks(int n_electrons,
    const integrals::LibcintInput& basis,
    const Eigen::Ref<const Eigen::MatrixXd>& overlap,
    const integrals::AoIntegralInput& ao,
    const RksOptions& options,
    const Eigen::MatrixXd& initial_density) {
  const auto grid = make_grid(basis, options.grid);
  const LdaFunctional functional(options.functional);
  const LdaIntegrator integrator(basis, grid, functional);
  RksResult result;
  result.grid_points = static_cast<std::size_t>(grid.weights.size());
  bool first = true;
  const auto evaluate = [&](const Eigen::Ref<const Eigen::MatrixXd>& density) {
    auto xc = integrator.evaluate(density);
    Eigen::MatrixXd coulomb = integrals::contract_jk(density, ao);
    scf::Evaluation physical;
    physical.electronic_energy =
        (density.array() * ao.ao_core_hamiltonian_matrix.array()).sum() +
        0.5 * (density.array() * coulomb.array()).sum() + xc.energy;
    physical.fock = ao.ao_core_hamiltonian_matrix + coulomb + xc.potential;
    if (first) {
      result.initial_electronic_energy = physical.electronic_energy;
      first = false;
    }
    result.xc_energy = xc.energy;
    result.integrated_electrons = xc.electron_count;
    return physical;
  };
  result.state = scf::solve_restricted(n_electrons, overlap,
      ao.ao_core_hamiltonian_matrix, evaluate, options.scf, initial_density);
  return result;
}

}  // namespace xmvb::dft
