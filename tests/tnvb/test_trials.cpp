#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include "input/loading/loader.hpp"
#include "tnvb/electronic/evaluator.hpp"
#include "vbscf/optimization/objective/function.hpp"

namespace {
void require(bool value,const char* message) {
  if (!value) throw std::runtime_error(message);
}
}

int main(int argc,char** argv) {
  try {
    if (argc!=2) throw std::invalid_argument("expected benzene fixture");
    auto loaded=xmvb::vb::load_vbscf_input_with_timings(argv[1]);
    xmvb::vb::SparseParameterLayout layout(loaded.input.orbital_preparation_input);
    const auto initial=layout.pack(loaded.input.orbital_preparation_input);
    xmvb::vb::OrbitalGradientEvaluator pullback;
    auto backend=std::make_shared<xmvb::tnvb::ElectronicEvaluator>();
    xmvb::vb::VbScfObjective objective(std::move(loaded.input),layout,{0},{1.0},
        loaded.nuclear_repulsion_energy,xmvb::vb::StructureEigensolver::Dense,{},&pullback,backend);
    Eigen::VectorXd gradient;
    const double accepted_energy=objective(initial,gradient);
    const auto accepted_orbitals=objective.input().orbital_preparation_input.orbital_value_table;
    const auto accepted_matrices=objective.gradient_result().scf_result.structure_matrices.hamiltonian_matrix;
    const auto verify_accepted=[&]() {
      require(objective.input().orbital_preparation_input.orbital_value_table==accepted_orbitals,"trial overwrote accepted orbitals");
      require(objective.gradient_result().scf_result.structure_matrices.hamiltonian_matrix==accepted_matrices,"trial overwrote accepted H/S result");
      require(objective.gradient_result().scf_result.total_energy==accepted_energy && objective.call_count()==1,"trial overwrote accepted energy/history");
      require(!objective.second_order_context() && objective.input().structure_data.alpha_det.empty(),"trial created determinant objects/context");
    };
    // Pair-member reversal leaves the C++ singlet CSF unchanged. It must not
    // silently change a structure's H/S row phase or its reported coefficient.
    auto reversed=objective.input();
    const int first_active=reversed.raw_structure_data.n_total_electrons-
        reversed.raw_structure_data.n_active_electrons;
    auto& raw=reversed.raw_structure_data.raw_structure_orbitals;
    std::swap(raw.at(first_active),raw.at(first_active+1));
    xmvb::vb::ElectronicStateRequest request;
    request.nuclear_repulsion_energy=loaded.nuclear_repulsion_energy;
    const auto reversed_result=backend->evaluate_energy_result(reversed,request).scf_result;
    const auto& original_result=objective.gradient_result().scf_result;
    const auto close=[](const auto& a,const auto& b) {
      if (a.size()!=b.size()) return false;
      for (std::size_t i=0;i<a.size();++i) if (std::abs(a[i]-b[i])>1e-10) return false;
      return true;
    };
    require(close(reversed_result.structure_matrices.hamiltonian_matrix,original_result.structure_matrices.hamiltonian_matrix) &&
        close(reversed_result.structure_matrices.overlap_matrix,original_result.structure_matrices.overlap_matrix),
        "singlet pair reversal changed the C++ structure phase");
    const int n=original_result.n_structures;
    double phase_overlap=0;
    for (int i=0;i<n;++i) phase_overlap+=reversed_result.eigenvector_matrix[i]*original_result.eigenvector_matrix[i];
    const double phase=phase_overlap<0?-1.0:1.0;
    for (int i=0;i<n;++i) require(std::abs(phase*reversed_result.eigenvector_matrix[i]-original_result.eigenvector_matrix[i])<1e-8,
        "singlet pair reversal changed structure coefficients");
    verify_accepted();
    Eigen::VectorXd trial=initial;
    trial[0]+=1e-4;
    // An ENERGY trial must be usable without a total MPS or any 3-PDM call.
    setenv("XMVB_BLOCK2_FORBID_3PDM","1",1);
    auto value=objective.evaluate_trial_energy(trial);
    require(value.valid && !value.gradient_ready && value.electronic_forward &&
        !value.forward_evaluation.has_value(),"ENERGY trial used determinant/gradient forward state");
    verify_accepted();
    // Fail inside the actual child worker, after orbital/integral preparation
    // and H/S construction. This is different from rejecting NaN up front:
    // neither a failed worker nor its trial cache may commit accepted results.
    bool worker_failed=false;
    try { (void)objective.evaluate_trial(trial); }
    catch (const std::runtime_error& error) {
      worker_failed=std::string(error.what()).find("TNVB worker failed")!=std::string::npos;
    }
    require(worker_failed,"forced FULL worker failure was not propagated");
    verify_accepted();
    unsetenv("XMVB_BLOCK2_FORBID_3PDM");
    objective.complete_trial(&value);
    require(value.gradient_ready && !value.gradient_result.second_order_context,"lazy TNVB completion failed");
    verify_accepted();
    auto bad=initial;
    bad[0]=std::numeric_limits<double>::quiet_NaN();
    bool rejected=false;
    try { (void)objective.evaluate_trial(bad,true); } catch (const std::exception&) { rejected=true; }
    require(rejected,"nonfinite trial was accepted");
    verify_accepted();
    objective.commit(std::move(value));
    require(objective.call_count()==2 && objective.input().orbital_preparation_input.orbital_value_table!=accepted_orbitals,"explicit trial commit failed");
    std::cout << "TNVB pair-phase, worker-failure and trial/accepted-state isolation passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n'; return 1;
  }
}
