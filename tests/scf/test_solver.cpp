#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "integrals/ao/graph.hpp"
#include "integrals/ao/input.hpp"
#include "integrals/ao/ri.hpp"
#include "scf/driver/solver.hpp"
#include "scf/hf/potential.hpp"
#include "scf/mixing/diis.hpp"

namespace {
using Matrix = Eigen::MatrixXd;
using xmvb::integrals::AoIntegralInput;
using xmvb::integrals::RiAoFactorization;
using xmvb::scf::Evaluation;
using xmvb::scf::Options;
using xmvb::scf::Result;

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

bool throws_invalid_argument(const std::function<void()>& action) {
  try {
    action();
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

/** @brief Verify closed-shell density, orbitals, and final evaluation agree. */
void check_result(const Result& result, int ne, const Matrix& overlap,
                  const Matrix& core,
                  const xmvb::scf::Evaluate& evaluate,
                  double tolerance, const std::string& label) {
  require(result.converged, label + ": expected convergence");
  require(result.iterations > 0, label + ": iteration count must be positive");
  require(result.density.rows() == overlap.rows() &&
              result.fock.rows() == overlap.rows() &&
              result.orbitals.rows() == overlap.rows(),
          label + ": result dimensions are inconsistent");
  require(std::abs((result.density * overlap).trace() - ne) < tolerance,
          label + ": Tr(P S) does not equal electron count");
  require((result.density * overlap * result.density -
           2.0 * result.density).norm() < tolerance,
          label + ": P S P != 2 P for a restricted determinant");
  require((result.density - 2.0 * result.orbitals.leftCols(ne / 2) *
               result.orbitals.leftCols(ne / 2).transpose()).norm() < tolerance,
          label + ": density does not match occupied orbitals");

  const Evaluation final = evaluate(result.density);
  require((result.fock - final.fock).norm() < tolerance,
          label + ": returned Fock matrix is not evaluated at returned density");
  require(std::abs(result.electronic_energy - final.electronic_energy) < tolerance,
          label + ": returned energy is not evaluated at returned density");
  require((result.fock - result.fock.transpose()).norm() < tolerance,
          label + ": returned Fock matrix is not symmetric");
  require((result.fock * result.density * overlap -
           overlap * result.density * result.fock).norm() < 10.0 * tolerance,
          label + ": final density has a nonzero F P S - S P F residual");

  // The generalized eigenproblem is checked directly, without an assumed
  // canonical orbital phase or any system-specific reference energies.
  const Matrix residual = result.fock * result.orbitals -
      overlap * result.orbitals * result.orbital_energies.asDiagonal();
  require(residual.norm() < 10.0 * tolerance,
          label + ": returned orbitals do not diagonalize the returned Fock matrix");
  (void)core;
}

Matrix make_overlap() {
  Matrix overlap(3, 3);
  overlap << 1.0, 0.18, -0.06,
             0.18, 1.25, 0.11,
             -0.06, 0.11, 0.92;
  return overlap;
}

Matrix make_core() {
  Matrix core(3, 3);
  core << -1.1, 0.07, -0.03,
          0.07, -0.45, 0.05,
          -0.03, 0.05, 0.16;
  return core;
}

RiAoFactorization make_ri() {
  constexpr int n = 3;
  constexpr int npair = n * (n + 1) / 2;
  RiAoFactorization ri;
  ri.n_basis_functions = n;
  ri.n_auxiliary_functions = 4;
  ri.n_packed_ao_pairs = npair;
  ri.metric_whitened_ao_pair_factors.resize(4, npair);
  for (int a = 0; a < 4; ++a) for (int p = 0; p < npair; ++p)
    ri.metric_whitened_ao_pair_factors(a, p) =
        0.12 * std::sin(0.31 * (a + 1) * (p + 2)) +
        0.035 * (a + 1) * (p + 1);
  return ri;
}

AoIntegralInput make_exact_from_ri(const RiAoFactorization& ri) {
  const int n = ri.n_basis_functions;
  const int npair = ri.n_packed_ao_pairs;
  std::vector<int> left, right;
  std::vector<double> values;
  const Matrix pair_eri = ri.metric_whitened_ao_pair_factors.transpose() *
                          ri.metric_whitened_ao_pair_factors;
  for (int p = 0; p < npair; ++p) for (int q = 0; q <= p; ++q) {
    left.push_back(p);
    right.push_back(q);
    values.push_back(pair_eri(p, q));
  }
  AoIntegralInput ao;
  ao.n_basis_functions = n;
  ao.pair_graph = xmvb::integrals::build_ao_pair_graph(
      std::move(left), std::move(right), std::move(values), n);
  return ao;
}

xmvb::scf::Evaluate make_hf_model(const Matrix& core,
                                  const AoIntegralInput& ao) {
  return [&core, &ao](const Eigen::Ref<const Matrix>& density) {
    (void)core;
    return xmvb::scf::evaluate_hf(density, ao);
  };
}

xmvb::scf::Evaluate make_hf_model(const Matrix& core,
                                  const RiAoFactorization& ri) {
  return [&core, &ri](const Eigen::Ref<const Matrix>& density) {
    return xmvb::scf::evaluate_hf(density, core, ri);
  };
}

void check_constant_model(const Matrix& overlap, const Matrix& core) {
  const xmvb::scf::Evaluate constant = [&core](
      const Eigen::Ref<const Matrix>& density) {
    Evaluation evaluation;
    evaluation.fock = core;
    evaluation.electronic_energy = (density.array() * core.array()).sum();
    return evaluation;
  };
  Options options;
  options.density_tolerance = 1e-11;
  const Result result = xmvb::scf::solve_restricted(
      2, overlap, core, constant, options);
  check_result(result, 2, overlap, core, constant, 1e-9, "constant model");
}

void check_nonlinear_models(const Matrix& overlap, const Matrix& core) {
  const RiAoFactorization ri = make_ri();
  AoIntegralInput ao = make_exact_from_ri(ri);
  ao.ao_core_hamiltonian_matrix = core;
  Matrix sample_density(3, 3);
  sample_density << 0.8, -0.1, 0.04,
                    -0.1, 0.3, 0.09,
                    0.04, 0.09, 0.6;
  const Evaluation exact_sample = xmvb::scf::evaluate_hf(sample_density, ao);
  const Evaluation ri_sample = xmvb::scf::evaluate_hf(sample_density, core, ri);
  require((exact_sample.fock - ri_sample.fock).norm() < 1e-12 &&
              std::abs(exact_sample.electronic_energy -
                       ri_sample.electronic_energy) < 1e-12,
          "exact and RI HF evaluators disagree for shared synthetic ERIs");
  const auto exact_model = make_hf_model(core, ao);
  const auto ri_model = make_hf_model(core, ri);
  Matrix direction(3, 3);
  direction << 0.12, -0.17, 0.09,
              -0.17, 0.04, 0.11,
               0.09, 0.11, -0.08;
  constexpr double step = 1e-5;
  for (const auto& model : {exact_model, ri_model}) {
    const double energy_derivative =
        (model(sample_density + step * direction).electronic_energy -
         model(sample_density - step * direction).electronic_energy) / (2 * step);
    require(std::abs(energy_derivative -
        (model(sample_density).fock.array() * direction.array()).sum()) < 1e-9,
        "HF energy derivative does not equal the physical Fock contraction");
  }
  Options options;
  options.density_tolerance = 1e-9;
  const Result exact = xmvb::scf::solve_restricted(
      2, overlap, core, exact_model, options);
  const Result fitted = xmvb::scf::solve_restricted(
      2, overlap, core, ri_model, options);
  check_result(exact, 2, overlap, core, exact_model, 2e-8,
               "nonlinear exact-ERI model");
  check_result(fitted, 2, overlap, core, ri_model, 2e-8,
               "nonlinear RI model");
  require((exact.density - fitted.density).norm() < 2e-8,
          "exact and RI models built from shared factors converged differently");
  require(std::abs(exact.electronic_energy - fitted.electronic_energy) < 2e-8,
          "exact and RI energies from shared factors differ");

  // One iteration cannot converge from a non-self-consistent starting density.
  Matrix initial(3, 3);
  initial << 0.4, 0.03, -0.01,
             0.03, 0.8, 0.02,
             -0.01, 0.02, 0.8;
  Options bounded = options;
  bounded.max_iterations = 1;
  bounded.density_tolerance = 1e-14;
  const Result limited = xmvb::scf::solve_restricted(
      2, overlap, core, exact_model, bounded, initial);
  require(!limited.converged,
          "iteration-limited SCF incorrectly reported convergence");
  require(limited.iterations == 1,
          "iteration-limited SCF reported the wrong iteration count");

  bounded.max_iterations = 0;
  require(throws_invalid_argument([&] {
    xmvb::scf::solve_restricted(
        2, overlap, core, exact_model, bounded, initial);
  }), "zero iteration limit was not rejected");
}

void check_invalid_inputs(const Matrix& overlap, const Matrix& core) {
  const xmvb::scf::Evaluate evaluator = [&core](
      const Eigen::Ref<const Matrix>& density) {
    return Evaluation{core, (density.array() * core.array()).sum()};
  };
  require(throws_invalid_argument([&] {
    xmvb::scf::solve_restricted(1, overlap, core, evaluator);
  }), "odd restricted electron count was not rejected");
  require(throws_invalid_argument([&] {
    xmvb::scf::solve_restricted(8, overlap, core, evaluator);
  }), "electron count above basis capacity was not rejected");
  Matrix wrong_size = Matrix::Identity(2, 2);
  require(throws_invalid_argument([&] {
    xmvb::scf::solve_restricted(2, wrong_size, core, evaluator);
  }), "mismatched overlap/core dimensions were not rejected");
}

/** @brief DIIS depends on residual directions, not their absolute units. */
void check_diis_scaling() {
  Matrix f0 = Matrix::Zero(2, 2);
  Matrix f1 = Matrix::Identity(2, 2);
  Matrix e0(2, 2), e1(2, 2);
  e0 << 1, 0, 0, 0;
  e1 << 0, 0, 0, 2;
  const Matrix expected = 0.2 * f1;
  for (double scale : {1.0, 1e-14, 1e-150, 1e150}) {
    xmvb::scf::Diis diis(6);
    diis.push(f0, scale * e0);
    diis.push(f1, scale * e1);
    Matrix mixed;
    require(diis.extrapolate(&mixed) && (mixed - expected).norm() < 1e-12,
            "DIIS changed under a common residual scaling");
    // Duplicating a residual leaves the constrained residual minimizer intact.
    diis.push(f1, scale * e1);
    require(diis.extrapolate(&mixed) && (mixed - expected).norm() < 1e-12,
            "dependent DIIS history disabled or changed extrapolation");
  }
  xmvb::scf::Diis zero(2);
  zero.push(f0, Matrix::Zero(2, 2));
  zero.push(f1, Matrix::Zero(2, 2));
  Matrix mixed;
  require(zero.extrapolate(&mixed) && (mixed - f1).norm() == 0.0,
          "zero-error DIIS history did not retain the latest Fock");
}

}  // namespace

int main() {
  try {
    const Matrix overlap = make_overlap();
    const Matrix core = make_core();
    check_constant_model(overlap, core);
    check_nonlinear_models(overlap, core);
    check_invalid_inputs(overlap, core);
    check_diis_scaling();
    std::cout << "restricted SCF invariants and exact/RI nonlinear models: PASS\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
