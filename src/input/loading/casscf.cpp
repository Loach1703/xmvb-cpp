#include "input/loading/casscf.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

extern "C" {
#include "cint.h"
}

#include "casscf/orbitals/frame.hpp"
#include "guess/block.hpp"
#include "guess/initial.hpp"
#include "input/deck/model.hpp"
#include "input/deck/primary_basis.hpp"
#include "integrals/ao/libcint.hpp"
#include "libcint/direct_shell.hpp"
#include "libcint/materialized_input.hpp"
#include "libcint/materialized_provider.hpp"
#include "vbscf/orbitals/preparation/input.hpp"

namespace xmvb::input {
namespace {

Eigen::MatrixXd build_overlap(const integrals::LibcintInput& libcint_input) {
  libcint::LibcintDirectShellEvaluator evaluator(libcint_input);
  const int n=evaluator.n_basis_functions();
  Eigen::MatrixXd overlap=Eigen::MatrixXd::Zero(n,n);
  for (int right=0; right<libcint_input.n_shells; ++right)
    for (int left=0; left<=right; ++left) {
      const auto block=evaluator.evaluate_overlap_shell_pair(left,right);
      for (int j=0; j<block.right_ao_count; ++j)
        for (int i=0; i<block.left_ao_count; ++i) {
          const int row=block.left_ao_offset+i;
          const int col=block.right_ao_offset+j;
          const double value=block.values[i+j*block.left_ao_count];
          overlap(row,col)=value;
          overlap(col,row)=value;
        }
    }
  return overlap;
}

std::pair<int,double> electron_count_and_nuclear_energy(
    const integrals::LibcintInput& libcint_input, int charge) {
  integrals::validate_libcint_input_shape(libcint_input);
  int electrons=-charge;
  double nuclear=0.0;
  for (int i=0; i<libcint_input.n_atoms; ++i) {
    const int z=libcint_input.atm[i*ATM_SLOTS+CHARGE_OF];
    const int offset=libcint_input.atm[i*ATM_SLOTS+PTR_COORD];
    electrons+=z;
    for (int j=0; j<i; ++j) {
      const int other=libcint_input.atm[j*ATM_SLOTS+PTR_COORD];
      const double dx=libcint_input.env[offset]-libcint_input.env[other];
      const double dy=libcint_input.env[offset+1]-libcint_input.env[other+1];
      const double dz=libcint_input.env[offset+2]-libcint_input.env[other+2];
      const double radius=std::sqrt(dx*dx+dy*dy+dz*dz);
      if (!(radius>0.0))
        throw std::invalid_argument("CASSCF geometry contains overlapping nuclei");
      nuclear+=static_cast<double>(z)*
          libcint_input.atm[j*ATM_SLOTS+CHARGE_OF]/radius;
    }
  }
  if (electrons<1)
    throw std::invalid_argument("CASSCF molecular electron count is not positive");
  return {electrons,nuclear};
}

Eigen::MatrixXd build_seed(const vb::InputDeck& deck,
    const std::string& input_path, const CasscfInput& input) {
  const int n=input.ao.n_basis_functions;
  const int occupied=input.ncore+input.nactive;
  vb::OrbitalPreparationInput bridge;
  bridge.n_basis_functions=n;
  bridge.n_orbitals=occupied;
  bridge.n_active_orbitals=input.nactive;
  bridge.n_total_electrons=input.n_electrons;
  bridge.n_active_electrons=input.n_electrons-2*input.ncore;
  bridge.spin_multiplicity=input.twice_spin+1;
  bridge.orbital_type=deck.metadata.orbital_type;
  bridge.ao_overlap_matrix=input.overlap;
  bridge.ao_normalization=vb::build_ao_normalization(input.overlap);
  bridge.orbital_basis_counts.assign(occupied,n);
  bridge.original_orbital_basis_counts.assign(occupied,n);
  bridge.orbital_basis_index_table.resize(static_cast<std::size_t>(occupied)*n);
  for (int orbital=0; orbital<occupied; ++orbital)
    for (int mu=0; mu<n; ++mu)
      bridge.orbital_basis_index_table[static_cast<std::size_t>(orbital)*n+mu]=mu+1;
  vb::build_initial_orbital_guess(
      input_path,deck.metadata.guess_type,input.libcint_input,input.ao,nullptr,
      deck.guess_block.present ? &deck.guess_block.raw_lines : nullptr,&bridge);
  Eigen::MatrixXd seed(n,occupied);
  for (int orbital=0; orbital<occupied; ++orbital)
    for (int mu=0; mu<n; ++mu)
      seed(mu,orbital)=bridge.orbital_value_table[
          static_cast<std::size_t>(orbital)*n+mu];
  return seed;
}

}  // namespace

CasscfInput load_casscf_input(const vb::InputDeck& deck,
    const std::string& input_path) {
  const auto& metadata=deck.metadata;
  if (metadata.request_ri_two_electron_mode)
    throw std::invalid_argument("CASSCF currently requires exact AO integrals; RI requested");
  if (metadata.state_average_count!=1)
    throw std::invalid_argument("CASSCF currently supports one state only");
  if (metadata.declared_active_orbitals<1 ||
      metadata.declared_active_electrons<1 ||
      metadata.declared_spin_multiplicity<1)
    throw std::invalid_argument("CASSCF requires NAO, NAE, and NMUL in $CTRL");
  const auto basis=vb::build_input_deck_basis(deck,metadata.basis_name);
  CasscfInput result;
  result.libcint_input=basis.libcint_input;
  result.basis=basis.basis_display_name;
  result.overlap=build_overlap(result.libcint_input);
  const auto [total_electrons,nuclear]=
      electron_count_and_nuclear_energy(result.libcint_input,metadata.total_charge);
  result.nactive=metadata.declared_active_orbitals;
  result.n_electrons=total_electrons;
  result.twice_spin=metadata.declared_spin_multiplicity-1;
  result.nuclear_energy=nuclear;
  result.max_iterations=metadata.requested_scf_max_iterations>0 ?
      metadata.requested_scf_max_iterations : 2000;
  result.request_molden_output=metadata.request_molden_output;
  const int inactive_electrons=total_electrons-metadata.declared_active_electrons;
  if (inactive_electrons<0 || inactive_electrons%2!=0)
    throw std::invalid_argument("CASSCF NAE leaves an invalid closed-shell inactive count");
  result.ncore=inactive_electrons/2;
  const int active_electrons=metadata.declared_active_electrons;
  if (result.ncore+result.nactive>result.overlap.rows() ||
      active_electrons>2*result.nactive ||
      result.twice_spin>active_electrons ||
      result.twice_spin>2*result.nactive-active_electrons ||
      (active_electrons-result.twice_spin)%2!=0)
    throw std::invalid_argument("CASSCF active electron, spin, or orbital partition is invalid");
  libcint::LibcintMaterializedIntegralProvider provider;
  result.ao=libcint::build_materialized_ao_integral_input(
      provider.build(result.libcint_input));
  const Eigen::MatrixXd seed=build_seed(deck,input_path,result);
  result.orbitals=casscf::complete_metric_orthonormal_frame(
      seed,result.overlap,result.ncore,result.nactive);
  return result;
}

}  // namespace xmvb::input
