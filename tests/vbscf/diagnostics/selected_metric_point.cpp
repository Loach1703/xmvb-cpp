#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#include "input/loading/loader.hpp"
#include "linalg/metric_solve.hpp"
#include "linalg/eigensolver.hpp"
#include "output/trace/binary.hpp"
#include "vbscf/derivatives/gradient/orbital/evaluator.hpp"
#include "vbscf/derivatives/hessian/context/accepted_point.hpp"
#include "vbscf/orbitals/preparation/cholesky.hpp"
#include "vbscf/structures/overlap/representation.hpp"
#include "vbscf/structures/overlap/spin_metric.hpp"
#include "vbscf/derivatives/hessian/coupled/structure.hpp"

namespace {

using xmvb::vb::RepresentationOverlap;

/** @brief Same current H/M action, seeds and targets; no orbital optimization. */
void inspect_ci_models(const xmvb::vb::VbScfInput& input,
    const std::shared_ptr<xmvb::vb::AcceptedPointContext>& context) {
  const auto& action=*context->structure_action;
  const auto& spin=context->spin_metric;
  if (!spin) throw std::runtime_error("CI model inspection needs selected spin factors");
  const auto& diagonal=action.preconditioner_diagonal();
  const double energy=context->selected_state_energies.at(0);
  const Eigen::VectorXd denominator=diagonal.hamiltonian-energy*diagonal.overlap;
  std::cout << "model_dimension=" << action.n_structures() << '\n'
            << "model_negative_diagonals=" << (denominator.array()<0.0).count() << '\n'
            << "model_denominator_min=" << denominator.minCoeff() << '\n'
            << "model_denominator_max=" << denominator.maxCoeff() << '\n';
  // Reconstruct the discarded trace model only in this diagnostic. Group
  // original selected columns by occupation; never form a physical matrix.
  std::map<std::vector<int>,std::vector<int>> groups;
  std::vector<bool> seen(action.n_structures(),false);
  const auto& expansion=input.structure_data;
  for (std::size_t p=0;p<expansion.determinant_to_structure_terms.size();++p) {
    std::vector<int> occupation=expansion.alpha_det[p];
    occupation.insert(occupation.end(),expansion.beta_det[p].begin(),expansion.beta_det[p].end());
    std::sort(occupation.begin(),occupation.end());
    for (const auto& term:expansion.determinant_to_structure_terms[p])
      if (!seen[term.structure_index] && term.coefficient!=0.0) {
        seen[term.structure_index]=true;
        groups[occupation].push_back(term.structure_index);
      }
  }
  Eigen::VectorXd trace_inverse(action.n_structures());
  int negative_blocks=0;
  double worst_spread=1.0;
  for (const auto& [occupation,indices]:groups) {
    double g=0.0,h=0.0,m=0.0,lo=std::numeric_limits<double>::infinity(),hi=0.0;
    for (int i:indices) {
      g+=spin->diagonal()[i]; h+=diagonal.hamiltonian[i]; m+=diagonal.overlap[i];
      const double value=std::abs(denominator[i])/spin->diagonal()[i];
      lo=std::min(lo,value); hi=std::max(hi,value);
    }
    if (h-energy*m<0.0) ++negative_blocks;
    if (lo>0.0) worst_spread=std::max(worst_spread,hi/lo);
    for (int i:indices) trace_inverse[i]=g; // Energy-dependent scale filled below.
  }
  std::cout << "trace_model_negative_blocks=" << negative_blocks << '\n'
            << "trace_model_worst_within_block_spread=" << worst_spread << '\n';
  const auto accuracy=context->structure_solve_accuracy;
  const xmvb::linalg::GeneralizedEigenAction images=[&](const Eigen::Ref<const Eigen::MatrixXd>& x) {
    auto pair=action.apply(x);
    return xmvb::linalg::GeneralizedEigenActionResult{
        std::move(pair.hamiltonian),std::move(pair.overlap)};
  };
  for (int model:{0,1,2}) {
    auto options=xmvb::linalg::make_davidson_options(action.n_structures(),1,
        accuracy.energy_tolerance,accuracy.gradient_tolerance,accuracy.stationarity_tolerance());
    options.residual_preconditioner=[&](const Eigen::VectorXd& r,double e) {
      if (model==0) return spin->precondition(r,diagonal.overlap);
      if (model==1) {
        Eigen::VectorXd z=spin->solve(r).col(0);
        for (const auto& [occupation,indices]:groups) {
          double h=0.0,m=0.0;
          for (int i:indices) { h+=diagonal.hamiltonian[i]; m+=diagonal.overlap[i]; }
          double divisor=h-e*m;
          const double floor=std::sqrt(std::numeric_limits<double>::epsilon())*
              std::max({trace_inverse[indices[0]],std::abs(h),std::abs(e*m)});
          if (std::abs(divisor)<floor) divisor=std::copysign(floor,divisor==0.0 ? 1.0:divisor);
          for (int i:indices) z[i]*=trace_inverse[i]/divisor;
        }
        return z;
      }
      const Eigen::VectorXd floor=std::sqrt(std::numeric_limits<double>::epsilon())*
          diagonal.hamiltonian.cwiseAbs().cwiseMax((e*diagonal.overlap).cwiseAbs()).cwiseMax(1.0);
      const Eigen::VectorXd d=(diagonal.hamiltonian-e*diagonal.overlap).cwiseAbs().cwiseMax(floor);
      return Eigen::VectorXd(spin->solve_diagonal_model(r,d).col(0));
    };
    const auto start=std::chrono::steady_clock::now();
    const auto solved=xmvb::linalg::GeneralizedEigensolver{}.solve_davidson(
        images,diagonal.hamiltonian,diagonal.overlap,options);
    std::cout << "ci_model=" << (model==0 ? "metric_spin":model==1 ? "discarded_trace":"diagonal_spin")
              << ",actions=" << solved.block_actions
              << ",energy=" << solved.eigenpairs.eigenvalues[0]
              << ",residual=" << solved.stationarity_residual_norms[0]
              << ",wall=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
              << '\n';
  }
  const xmvb::vb::StructureTangentOperator tangent(context,action);
  const int size=tangent.tangent_size();
  Eigen::VectorXd x(size);
  for (int i=0;i<size;++i) x[i]=std::sin(0.37*(i+1))+std::cos(0.11*(i+1));
  x.normalize();
  const Eigen::VectorXd rhs=tangent.apply_hessian_coordinates(x);
  for (bool matched:{false,true}) {
    Eigen::VectorXd z;
    if (matched) z=tangent.apply_inverse_augmented_hessian_diagonal(rhs,0.0);
    else {
      Eigen::MatrixXd raw=tangent.expand(rhs).scaled_coefficients;
      Eigen::VectorXd d=denominator;
      for (Eigen::Index i=0;i<d.size();++i) if (std::abs(d[i])<1e-8) d[i]=1e-8;
      raw.col(0).array()/=d.array();
      z=tangent.project_coordinates(raw);
    }
    const Eigen::VectorXd az=tangent.apply_hessian_coordinates(z);
    const double scale=rhs.dot(az)/az.squaredNorm();
    std::cout << "ah_model=" << (matched ? "diagonal_spin":"raw_diagonal")
              << ",inverse_work=" << rhs.dot(z)
              << ",best_linear_residual=" << (rhs-scale*az).norm()/rhs.norm() << '\n';
  }
}

void print_spectrum(const char* label, const Eigen::MatrixXd& matrix) {
  const Eigen::MatrixXd symmetric = 0.5 * (matrix + matrix.transpose()).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(
      symmetric, Eigen::EigenvaluesOnly);
  if (eigensolver.info() != Eigen::Success)
    throw std::runtime_error(std::string(label) + " eigensolver failed");
  const double smallest = eigensolver.eigenvalues()[0];
  const double largest = eigensolver.eigenvalues().tail(1)[0];
  std::cout << label << "_lambda_min=" << smallest << '\n'
            << label << "_lambda_max=" << largest << '\n'
            << label << "_condition="
            << (smallest > 0.0 ? largest / smallest
                               : std::numeric_limits<double>::infinity()) << '\n';
}

template<class Action>
void inspect_rhs(const char* label, const Eigen::VectorXd& rhs,
                 const Action& action, const Eigen::VectorXd& jacobi,
                 const Eigen::LLT<Eigen::MatrixXd>& factor) {
  const double rhs_norm = rhs.norm();
  const double target = std::sqrt(std::numeric_limits<double>::epsilon()) * rhs_norm;
  std::cout << label << "_rhs_norm=" << rhs_norm << '\n'
            << label << "_target=" << target << '\n';

  xmvb::linalg::MetricSolveStats stats;
  try {
    const Eigen::VectorXd solution = xmvb::linalg::solve_metric(
        rhs, action,
        [&](const Eigen::VectorXd& residual) {
          return Eigen::VectorXd(residual.cwiseQuotient(jacobi));
        }, &stats);
    const double true_residual = (rhs - action(solution)).norm();
    std::cout << label << "_cg_status=certified\n"
              << label << "_cg_solution_norm=" << solution.norm() << '\n'
              << label << "_cg_true_residual=" << true_residual << '\n'
              << label << "_cg_relative_residual="
              << (rhs_norm > 0.0 ? true_residual / rhs_norm : 0.0) << '\n';
  } catch (const std::exception& error) {
    std::cout << label << "_cg_status=failed\n"
              << label << "_cg_failure=" << error.what() << '\n';
  }
  std::cout << label << "_cg_iterations=" << stats.iterations << '\n'
            << label << "_cg_actions=" << stats.actions << '\n'
            << label << "_cg_certifications=" << stats.certifications << '\n'
            << label << "_cg_last_certified_relative_residual="
            << stats.relative_residual << '\n';

  if (factor.info() != Eigen::Success) {
    std::cout << label << "_llt_status=failed\n";
    return;
  }
  const Eigen::VectorXd solution = factor.solve(rhs);
  if (factor.info() != Eigen::Success || !solution.allFinite()) {
    std::cout << label << "_llt_status=failed\n";
    return;
  }
  const double true_residual = (rhs - action(solution)).norm();
  std::cout << label << "_llt_status="
            << (true_residual <= target ? "certified" : "above_tolerance") << '\n'
            << label << "_llt_solution_norm=" << solution.norm() << '\n'
            << label << "_llt_true_residual=" << true_residual << '\n'
            << label << "_llt_relative_residual="
            << (rhs_norm > 0.0 ? true_residual / rhs_norm : 0.0) << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2 || argc > 3)
      throw std::invalid_argument(
          "usage: test_selected_metric_point INPUT.xmi [OUTPUT_PREFIX|--ci-model]");
    auto loaded = xmvb::vb::load_vbscf_input_with_timings(argv[1]);
    const Eigen::MatrixXd no_initial_eigenvectors;
    xmvb::vb::OrbitalGradientEvaluator evaluator;
    const auto evaluated = evaluator.evaluate_without_reference_energy_gradient(
        loaded.input, {0}, {1.0}, loaded.nuclear_repulsion_energy,
        xmvb::vb::StructureEigensolver::Davidson,
        xmvb::vb::StructureSolveAccuracy{}, no_initial_eigenvectors);
    const auto& context = evaluated.second_order_context;
    if (!context || !context->structure_action)
      throw std::runtime_error("accepted structure action is unavailable");
    const auto& action = *context->structure_action;
    if (argc==3 && std::string(argv[2])=="--ci-model") {
      std::cout << std::setprecision(17);
      inspect_ci_models(loaded.input,context);
      return 0;
    }
    const Eigen::MatrixXd& coefficients = context->selected_state_eigenvectors;
    const auto images = context->selected_state_structure_images
        ? *context->selected_state_structure_images : action.apply(coefficients);
    const auto& orbitals = loaded.input.orbital_preparation_input;
    const int n_active = orbitals.n_active_orbitals;
    const auto& overlap =
        context->prepared_active_space.orbital_result.active_orbital_overlap_matrix;
    const Eigen::MatrixXd upper = xmvb::vb::overlap_upper(
        Eigen::Map<const Eigen::MatrixXd>(overlap.data(), n_active, n_active));
    const auto& alpha =
        context->same_spin_pair_cache.alpha_reuse_table.unique_determinants;
    const auto& beta =
        context->same_spin_pair_cache.beta_reuse_table.unique_determinants;
    const auto prepared = RepresentationOverlap(action, alpha, beta, n_active)
                              .prepare_metric(upper);
    const int dimension = action.n_structures();
    std::cout << std::setprecision(17)
              << "direct_ci=" << action.supports_integral_direction() << '\n'
              << "dimension=" << dimension << '\n'
              << "resident_selected_matrix=" << (prepared.selected_matrix().size() != 0)
              << '\n';
    if (!RepresentationOverlap::PreparedMetric::fits_selected(dimension))
      throw std::runtime_error("selected Gram exceeds bounded diagnostic workspace");

    // The identity block is permitted only by the existing selected-matrix
    // budget. No full spin-carrier or determinant-product matrix is formed.
    const Eigen::MatrixXd metric = prepared.selected_matrix().size()
        ? prepared.selected_matrix()
        : prepared(Eigen::MatrixXd::Identity(dimension, dimension));
    const auto original_action = [&](const Eigen::VectorXd& vector) {
      return Eigen::VectorXd(prepared(vector).col(0));
    };
    const Eigen::MatrixXd symmetric = 0.5 * (metric + metric.transpose()).eval();
    std::cout << "metric_symmetry_relative="
              << (metric - metric.transpose()).norm() / std::max(1.0, metric.norm())
              << '\n';
    print_spectrum("metric", metric);
    const Eigen::VectorXd jacobi = action.preconditioner_diagonal().overlap;
    if (jacobi.size() != dimension || !jacobi.allFinite() ||
        (jacobi.array() <= 0.0).any())
      throw std::runtime_error("accepted Jacobi diagonal is invalid");
    const Eigen::VectorXd inverse_sqrt = jacobi.array().rsqrt();
    const Eigen::MatrixXd scaled =
        inverse_sqrt.asDiagonal() * symmetric * inverse_sqrt.asDiagonal();
    print_spectrum("jacobi_scaled", scaled);
    std::cout << "jacobi_metric_diagonal_relative_gap="
              << (jacobi - metric.diagonal()).norm() /
                     std::max(1.0, metric.diagonal().norm()) << '\n';

    const Eigen::MatrixXd state_h = coefficients.transpose() * images.hamiltonian;
    const Eigen::MatrixXd transport_rhs =
        2.0 * (images.hamiltonian - images.overlap * state_h);
    Eigen::VectorXd deterministic(dimension);
    for (int i = 0; i < dimension; ++i)
      deterministic[i] = std::sin(0.37 * (i + 1)) +
                         std::cos(0.11 * (i + 1));
    const Eigen::VectorXd deterministic_image = original_action(deterministic);
    std::cout << "deterministic_action_matrix_relative_gap="
              << (deterministic_image - metric * deterministic).norm() /
                     std::max(1.0, deterministic_image.norm()) << '\n';

    if (argc == 3) {
      const std::filesystem::path prefix(argv[2]);
      xmvb::output::write_binary_buffer(prefix.string() + ".metric.bin",
                                        metric.data(), metric.size());
      xmvb::output::write_binary_buffer(prefix.string() + ".gc0.bin",
                                        transport_rhs.data(), transport_rhs.size());
      xmvb::output::write_binary_buffer(prefix.string() + ".deterministic.bin",
                                        deterministic.data(), deterministic.size());
    }

    const Eigen::LLT<Eigen::MatrixXd> factor(symmetric);
    std::cout << "llt_factor_success=" << (factor.info() == Eigen::Success) << '\n';
    for (int state = 0; state < transport_rhs.cols(); ++state) {
      const std::string label = "transport_" + std::to_string(state);
      inspect_rhs(label.c_str(), transport_rhs.col(state), original_action,
                  jacobi, factor);
    }
    inspect_rhs("deterministic", deterministic, original_action, jacobi, factor);
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "selected_metric_point_failure=" << error.what() << '\n';
    return 1;
  }
}
