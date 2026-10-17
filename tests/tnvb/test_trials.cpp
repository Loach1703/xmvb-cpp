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
#include "vbscf/optimization/driver/optimizer.hpp"

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
    // A successful CLI exit at an iteration cap is not proof of a usable
    // optimizer: require an actual accepted P-LBFGS orbital step here.
    auto block_input=xmvb::vb::load_vbscf_input_with_timings(argv[1]);
    xmvb::vb::VbScfOptimizerOptions block_options;
    block_options.backend=xmvb::vb::VbScfOptimizerBackend::Plbfgs;
    block_options.structure_eigensolver=xmvb::vb::StructureEigensolver::Dense;
    block_options.max_iterations=1;
    block_options.verbose=false;
    block_options.retain_accepted_iteration_trace=true;
    const auto block=xmvb::vb::VbScfOptimizer{block_options,backend}.optimize(
        std::move(block_input.input),block_input.nuclear_repulsion_energy);
    require(block.n_iterations==1 && block.total_energy_history.size()==2,
        "TNVB block-LBFGS did not accept an orbital step");
    require(block.total_energy_history.back()<block.total_energy_history.front(),
        "TNVB block-LBFGS orbital step did not decrease energy");
    require(block.termination_reason=="max_iterations" || block.converged,
        "TNVB block-LBFGS stopped with an optimizer failure");
    require(block.accepted_iteration_trace.size()==2 &&
        block.accepted_iteration_trace.back().has_projected_gradient &&
        std::isfinite(block.final_projected_gradient_inf_norm),
        "TNVB block-LBFGS did not report its projected gradient");
    require(block.optimized_input.structure_data.alpha_det.empty() &&
        !block.optimized_input.complete_active_space && block.scf_result.n_structures==5,
        "TNVB block-LBFGS changed selected VB structures to determinants/CAS");
    require(block.scf_result.tnvb_diagnostics.has_value() &&
        block.scf_result.tnvb_diagnostics->fit_energy_error<1e-7,
        "TNVB block-LBFGS lost full-gradient tensor diagnostics");
    std::cout << "TNVB pair-phase, worker-failure and trial/accepted-state isolation passed\n";
    std::cout << "TNVB ISCF=8 accepted orbital step and selected-space invariants passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n'; return 1;
  }
}
