#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>

#include "vbscf/derivatives/gradient/orbital/evaluator.hpp"
#include "vbscf/derivatives/hessian/context/accepted_point.hpp"
#include "vbscf/integrals/ao/pairs/two_electron_index.hpp"
#include "vbscf/structures/expansion/expander.hpp"

namespace {

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

/** @brief Complete accepted-gradient producer, not only a StructureAction consumer. */
void benchmark_accepted_state16() {
  constexpr int n = 16;
  xmvb::vb::VbScfInput input;
  auto& orbitals = input.orbital_preparation_input;
  orbitals.n_basis_functions = orbitals.n_orbitals = orbitals.n_active_orbitals = n;
  orbitals.n_total_electrons = orbitals.n_active_electrons = n;
  orbitals.orbital_basis_counts.assign(n,n);
  orbitals.orbital_value_table.assign(n*n,0.0);
  orbitals.orbital_basis_index_table.resize(n*n);
  orbitals.ao_overlap_matrix = Eigen::MatrixXd::Identity(n,n);
  orbitals.ao_overlap_matrix.array() += 0.01;
  for (int i = 0; i < n; ++i) {
    orbitals.orbital_value_table[i*n+i] = 1.0;
    for (int j = 0; j < n; ++j) orbitals.orbital_basis_index_table[i*n+j] = j+1;
  }
  auto& ao = input.ao_integral_input;
  ao.n_basis_functions = n;
  ao.ao_core_hamiltonian_matrix = Eigen::MatrixXd::Zero(n,n);
  for (int i = 0; i < n; ++i) ao.ao_core_hamiltonian_matrix(i,i) = -0.9+0.03*i;
  std::vector<int> left,right;
  std::vector<double> values;
  for (int p = 0; p < n*(n+1)/2; ++p)
    for (int q = 0; q <= p; ++q) {
      left.push_back(p); right.push_back(q);
      values.push_back(0.02*std::cos(0.1*p)*std::cos(0.1*q)+
                       0.01*std::sin(0.07*p)*std::sin(0.07*q));
    }
  ao.pair_graph = xmvb::vb::build_ao_pair_graph(left,right,values,n);
  xmvb::vb::RawStructureData raw;
  raw.n_structures = 1;
  raw.n_total_electrons = raw.n_active_electrons = n;
  for (int i = 1; i <= n; ++i) raw.raw_structure_orbitals.push_back(i);
  input.structure_data = xmvb::vb::FullDeterminantStructureExpander().expand(raw);
  const auto start = std::chrono::steady_clock::now();
  const auto gradient = xmvb::vb::OrbitalGradientEvaluator().evaluate_without_reference_energy_gradient(
      input,{0},{1.0},0.0,xmvb::vb::StructureEigensolver::Davidson,{1e-7,1e-3},Eigen::MatrixXd{});
  const auto& accepted = *gradient.second_order_context;
  const auto& state = *accepted.structure_adjoint_state;
  require(accepted.selected_state_matrices.states.empty() && state.source_factors &&
              state.orthogonal_coefficients.empty() && state.residuals.empty(),
          "accepted-gradient producer retained a full product state");
  std::size_t raw_doubles = 0;
  for (const auto& s : accepted.selected_state_matrices.states) raw_doubles += s.coefficient_matrix.size();
  rusage usage{};
  getrusage(RUSAGE_SELF,&usage);
  const auto& cache = accepted.same_spin_pair_cache;
  const std::size_t removed_bytes = static_cast<std::size_t>(cache.n_unique_alpha)*cache.n_unique_beta*sizeof(double);
  require(std::isfinite(gradient.scf_result.total_energy) &&
              static_cast<std::size_t>(usage.ru_maxrss)*1024 < removed_bytes,
          "accepted-gradient peak RSS exceeds the removed raw-product allocation");
  std::cout << std::setprecision(17) << "orbitals=" << n << "\nstructures=1\n"
            << "alpha_strings=" << cache.n_unique_alpha << "\nbeta_strings=" << cache.n_unique_beta
            << "\naccepted_raw_product_doubles=" << raw_doubles << "\naccepted_orthogonal_product_doubles=0\naccepted_residual_product_doubles=0\n"
            << "removed_raw_product_bytes=" << removed_bytes << "\nenergy_Eh=" << gradient.scf_result.total_energy
            << "\ngradient_seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
            << "\nrss_mib=" << usage.ru_maxrss/1024.0 << '\n';
}

/** @brief Nonstationary fixed-state derivatives without any CI eigensolve. */
void check_keyframe(bool complete_carrier) {
  xmvb::vb::VbScfInput input;
  input.standard_two_electron_mode = xmvb::vb::StandardTwoElectronMode::Exact;
  auto& orbitals = input.orbital_preparation_input;
  orbitals.n_basis_functions = 3;
  orbitals.n_orbitals = orbitals.n_active_orbitals = 2;
  orbitals.n_total_electrons = orbitals.n_active_electrons = 2;
  orbitals.orbital_value_table = {1.0, 0.13, -0.17, 0.21, 0.94, 0.12};
  orbitals.orbital_basis_index_table = {1, 2, 3, 1, 2, 3};
  orbitals.orbital_basis_counts = {3, 3};
  orbitals.ao_overlap_matrix = Eigen::Matrix3d::Identity();
  orbitals.ao_overlap_matrix(0, 1) = orbitals.ao_overlap_matrix(1, 0) = 0.1;
  auto& ao = input.ao_integral_input;
  ao.n_basis_functions = 3;
  ao.ao_core_hamiltonian_matrix.resize(3, 3);
  ao.ao_core_hamiltonian_matrix << -1.1, 0.13, 0.04,
                                  0.13, -0.6, -0.08,
                                  0.04, -0.08, 0.2;
  std::vector<int> left, right;
  std::vector<double> values;
  for (int p = 0; p < 6; ++p) {
    for (int q = 0; q <= p; ++q) {
      left.push_back(p);
      right.push_back(q);
      values.push_back(0.07 * std::cos(0.2 * (p - q)));
    }
  }
  ao.pair_graph = xmvb::vb::build_ao_pair_graph(left, right, values, 3);
  auto& data = input.structure_data;
  for (int a = 0; a < (complete_carrier ? 2 : 1); ++a) {
    for (int b = 0; b < 2; ++b) {
      data.alpha_det.push_back({a});
      data.beta_det.push_back({b});
      data.determinant_to_structure_terms.push_back({{data.n_structures++, 1.0}});
    }
  }
  Eigen::MatrixXd c(data.n_structures, 2);
  for (int i = 0; i < c.rows(); ++i) {
    c(i, 0) = 0.37 + 0.19 * i;
    c(i, 1) = std::cos(0.73 * (i + 1));
  }
  const std::vector<int> states{0, data.n_structures - 1};
  const std::vector<double> weights{0.5, 0.5};
  const xmvb::vb::OrbitalGradientEvaluator evaluator;
  const auto relaxed = evaluator.evaluate(input, states, weights, 0.4);
  const auto& accepted = *relaxed.second_order_context;
  require(accepted.structure_action && accepted.structure_action->supports_integral_direction() &&
              accepted.selected_state_matrices.states.empty(),
          "direct-CI accepted point retained unused raw product states");
  const auto evaluate = [&](const Eigen::MatrixXd& coefficients) {
    return evaluator.evaluate_fixed_structure(input, states, weights, 0.4, coefficients);
  };
  const auto fixed = evaluate(c);
  require(fixed.gradient.scf_result.eigenvector_matrix.empty() &&
              fixed.gradient.scf_result.structure_overlap_diagonal.empty(),
          "fixed-state gradient constructed eigensolver output");
  const auto scaled = evaluate(2.7 * c);
  require(std::abs(scaled.gradient.scf_result.total_energy -
                   fixed.gradient.scf_result.total_energy) < 1e-11 &&
              (2.7 * scaled.coefficient_gradient - fixed.coefficient_gradient).norm() < 1e-11,
          "fixed-state normalization changed energy or residual");
  Eigen::Matrix2d change;
  change << 1.2, -0.3, 0.2, 0.8;
  const auto mixed = evaluate(c * change);
  require(std::abs(mixed.gradient.scf_result.total_energy -
                   fixed.gradient.scf_result.total_energy) < 1e-11 &&
              (mixed.coefficient_gradient * change.transpose() -
               fixed.coefficient_gradient).norm() < 1e-10,
          "equal-weight keyframe depends on the state frame");
  constexpr double delta = 1e-5;
  for (int i = 0; i < 6; ++i) {
    const double value = orbitals.orbital_value_table[i];
    orbitals.orbital_value_table[i] = value + delta;
    const double plus = evaluate(c).gradient.scf_result.total_energy;
    orbitals.orbital_value_table[i] = value - delta;
    const double minus = evaluate(c).gradient.scf_result.total_energy;
    orbitals.orbital_value_table[i] = value;
    require(std::abs((plus - minus) / (2 * delta) -
                     fixed.gradient.sparse_orbital_energy_gradient[i]) < 2e-7,
            "fixed-state orbital gradient differs from finite differences");
  }
  for (int state = 0; state < 2; ++state) {
    for (int i = 0; i < c.rows(); ++i) {
      const double value = c(i, state);
      c(i, state) = value + delta;
      const double plus = evaluate(c).gradient.scf_result.total_energy;
      c(i, state) = value - delta;
      const double minus = evaluate(c).gradient.scf_result.total_energy;
      c(i, state) = value;
      require(std::abs((plus - minus) / (2 * delta) -
                       fixed.coefficient_gradient(i, state)) < 2e-7,
              "fixed-state structure residual differs from finite differences");
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--accepted-state16") {
      benchmark_accepted_state16();
      return 0;
    }
    if (argc != 1) throw std::invalid_argument("usage: test_keyframe [--accepted-state16]");
    check_keyframe(true);
    check_keyframe(false);
    std::cout << "fixed-state keyframe actions and derivatives: passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
