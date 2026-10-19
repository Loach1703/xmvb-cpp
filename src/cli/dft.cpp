#include "cli/dft.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "cli/options.hpp"
#include "dft/driver/rks.hpp"
#include "dft/xc/functional.hpp"
#include "input/deck/model.hpp"
#include "input/loading/molecule.hpp"
#include "output/orbitals/writer.hpp"

namespace xmvb::cli {

int run_dft(const Options& command_line, const vb::InputDeck& deck) {
  const auto& metadata = deck.metadata;
  if (metadata.declared_spin_multiplicity != 1 || metadata.state_average_count != 1)
    throw std::invalid_argument("DFT currently supports closed-shell RKS, not open-shell or state averaging");
  if (metadata.request_ri_two_electron_mode ||
      command_line.load.standard_two_electron_mode == vb::StandardTwoElectronMode::ResolutionOfIdentity)
    throw std::invalid_argument("RKS currently requires exact AO integrals; RI will be added separately");
  if (metadata.guess_type != vb::kGuessTypeAuto || deck.guess_block.present)
    throw std::invalid_argument("RKS currently supports GUESS=AUTO (core Hamiltonian), not supplied guesses");
  if (command_line.optimizer_backend_explicit || command_line.structure_eigensolver_explicit ||
      metadata.scf_optimizer != vb::InputScfOptimizer::Unspecified ||
      !command_line.trace_directory.empty() || !command_line.ciah_trace_path.empty() ||
      !command_line.tnhvp_trace_path.empty() || !command_line.final_orbitals_path.empty() ||
      metadata.request_molden_output)
    throw std::invalid_argument("RKS uses the common SCF/DIIS driver; VB/CASSCF optimizer, trace and Molden options are unsupported");

  const auto start = std::chrono::steady_clock::now();
  auto molecule = input::load_molecule(deck);
  dft::RksOptions options;
  options.functional = metadata.xc_functional;
  options.grid.radial_points = metadata.grid_radial_points;
  options.grid.angular_points = metadata.grid_angular_points;
  options.scf.density_tolerance = metadata.scf_tolerance;
  if (metadata.requested_scf_max_iterations > 0)
    options.scf.max_iterations = metadata.requested_scf_max_iterations;
  if (command_line.max_iterations_explicit)
    options.scf.max_iterations = command_line.optimizer.max_iterations;

  std::cout << "Restricted Kohn-Sham DFT (LDA, exact J, SCF/DIIS)\n"
            << "Input: " << command_line.input_path << '\n'
            << "Basis: " << molecule.basis << '\n'
            << "AO functions: " << molecule.ao.n_basis_functions << '\n'
            << "Electrons: " << molecule.n_electrons << '\n'
            << "XC: " << options.functional << '\n'
            << "Libxc: " << dft::LdaFunctional::library_version() << '\n'
            << "Grid: " << options.grid.radial_points << " radial, "
            << options.grid.angular_points << " angular per atom\n"
            << "SCF density and residual tolerance: " << options.scf.density_tolerance << '\n';
  const auto result = dft::solve_rks(molecule.n_electrons, molecule.libcint_input,
      molecule.overlap, molecule.ao, options);
  if (result.state.converged) {
    const auto path = output::write_orbital_file(command_line.input_path, result.state.orbitals);
    std::cout << "Converged orbitals: " << path.string() << '\n';
  }
  std::cout << std::setprecision(15)
            << "Initial energy: " << result.initial_electronic_energy + molecule.nuclear_energy << " Eh\n"
            << "Final energy: " << result.state.electronic_energy + molecule.nuclear_energy << " Eh\n"
            << "XC energy: " << result.xc_energy << " Eh\n"
            << "Grid points: " << result.grid_points << '\n'
            << "Integrated electrons: " << result.integrated_electrons << '\n'
            << "SCF iterations: " << result.state.iterations << '\n'
            << "Converged: " << (result.state.converged ? "yes" : "no") << '\n'
            << "Total time: " << std::chrono::duration<double>(
                 std::chrono::steady_clock::now() - start).count() << " s\n";
  return result.state.converged ? 0 : 2;
}

}  // namespace xmvb::cli
