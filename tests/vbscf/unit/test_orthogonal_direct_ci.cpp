#include <array>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"
#include "vbscf/structures/orthogonal_ci/integrals.hpp"
#include "vbscf/structures/orthogonal_ci/metric.hpp"
#include "vbscf/structures/orthogonal_ci/planner.hpp"
#include "vbscf/structures/orthogonal_ci/sigma.hpp"
#include "vbscf/determinants/pairs/evaluator.hpp"
#include "vbscf/determinants/pairs/same_spin_cache.hpp"
#include "vbscf/integrals/active/two_electron/construction/indexer.hpp"
#include "vbscf/integrals/active/two_electron/construction/kernel.hpp"
#include "vbscf/structures/assembly/action.hpp"
#include "vbscf/structures/assembly/hamiltonian_overlap.hpp"

namespace {

void require(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void check_fused_response(
    const xmvb::vb::DirectCiSigmaAction& action,
    const Eigen::Ref<const Eigen::MatrixXd>& response,
    const Eigen::Ref<const Eigen::MatrixXd>& state) {
  const auto fused = action.apply_with_adjoint(response, state);
  const auto sigma = action.apply(response);
  const auto adjoint = action.integral_adjoint(response, state);
  require((fused.sigma - sigma).norm() < 3e-12 * std::max(1.0, sigma.norm()),
          "fused response sigma differs from separate action");
  require((fused.integrals.one_electron - adjoint.one_electron).norm() <
              3e-12 * std::max(1.0, adjoint.one_electron.norm()),
          "fused response one-electron adjoint differs from separate contraction");
  require((fused.integrals.pair_kernel - adjoint.pair_kernel).norm() <
              3e-12 * std::max(1.0, adjoint.pair_kernel.norm()),
          "fused response two-electron adjoint differs from separate contraction");
}

/** @brief Joint sweep agrees with independent H/S and differentiated adjoints. */
void check_joint_structure_response(const xmvb::vb::StructureAction& action,
    int n_alpha, int n_beta) {
  const int n = action.n_structures();
  xmvb::vb::SelectedStateDeterminantMatrices states;
  states.n_structures = n;
  states.n_unique_alpha = n_alpha;
  states.n_unique_beta = n_beta;
  states.normalized_state_weights = {0.3, 0.7};
  states.states.resize(2);
  Eigen::MatrixXd raw(n, 2);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < n; ++j) raw(j,i) = std::sin(0.23*(i+1)*(j+2));
    states.states[i].normalized_state_weight = states.normalized_state_weights[i];
    states.states[i].coefficient_matrix = action.expand_structure_block(raw.col(i));
  }
  const auto accepted = action.prepare_active_adjoint(states, {-1.3, -0.9});
  const Eigen::MatrixXd multipliers =
      (Eigen::Matrix2d() << 0.13, -0.04, 0.08, -0.21).finished();
  for (bool zero : {false, true}) {
    if (zero) raw.setZero();
    else raw = raw.array().cos().matrix();
    auto response = states;
    for (int i = 0; i < 2; ++i)
      response.states[i].coefficient_matrix = action.expand_structure_block(raw.col(i));
    auto sweep = action.apply_response(accepted, raw);
    const auto images = action.apply(raw);
    require((sweep.images.hamiltonian-images.hamiltonian).norm() <
        3e-11*std::max(1.0,images.hamiltonian.norm()) &&
        (sweep.images.overlap-images.overlap).norm() <
        3e-11*std::max(1.0,images.overlap.norm()), "joint response changed H/S images");
    const auto joint = action.finish_response(accepted, std::move(sweep), multipliers);
    const auto reference = action.active_integral_adjoint_direction(accepted,
        &response, &multipliers, {}, {}, {}, nullptr, false);
    require((joint.overlap-reference.overlap).norm() <
        3e-11*std::max(1.0,reference.overlap.norm()) &&
        (joint.one_electron-reference.one_electron).norm() <
        3e-11*std::max(1.0,reference.one_electron.norm()) &&
        (joint.pair_kernel-reference.pair_kernel).norm() <
        3e-11*std::max(1.0,reference.pair_kernel.norm()),
        "joint response changed the independent differentiated adjoint");
  }
}

void append_combinations(
    int n_orbitals,
    int n_electrons,
    int first,
    std::vector<int>* occupied,
    std::vector<std::vector<int>>* determinants) {
  if (static_cast<int>(occupied->size()) == n_electrons) {
    determinants->push_back(*occupied);
    return;
  }
  const int remaining = n_electrons - static_cast<int>(occupied->size());
  for (int orbital = first;
       orbital <= n_orbitals - remaining;
       ++orbital) {
    occupied->push_back(orbital);
    append_combinations(
        n_orbitals,
        n_electrons,
        orbital + 1,
        occupied,
        determinants);
    occupied->pop_back();
  }
}

std::vector<std::vector<int>> complete_space(
    int n_orbitals,
    int n_electrons) {
  std::vector<std::vector<int>> determinants;
  std::vector<int> occupied;
  append_combinations(
      n_orbitals,
      n_electrons,
      0,
      &occupied,
      &determinants);
  return determinants;
}

Eigen::MatrixXd explicit_exterior(
    const std::vector<std::vector<int>>& determinants,
    const Eigen::MatrixXd& orbital_transform) {
  const int dimension = static_cast<int>(determinants.size());
  const int n_electrons = static_cast<int>(determinants.front().size());
  Eigen::MatrixXd exterior(dimension, dimension);
  for (int target = 0; target < dimension; ++target) {
    for (int source = 0; source < dimension; ++source) {
      Eigen::MatrixXd minor(n_electrons, n_electrons);
      for (int row = 0; row < n_electrons; ++row) {
        for (int column = 0; column < n_electrons; ++column) {
          minor(row, column) = orbital_transform(
              determinants[target][row],
              determinants[source][column]);
        }
      }
      exterior(target, source) = minor.determinant();
    }
  }
  return exterior;
}

void check_exterior_transform(int width) {
  const auto determinants = complete_space(6, 3);
  Eigen::MatrixXd orbital = Eigen::MatrixXd::Zero(6, 6);
  for (int row = 0; row < 6; ++row) {
    orbital(row, row) = 0.8 + 0.07 * static_cast<double>(row + 1);
    for (int column = row + 1; column < 6; ++column) {
      orbital(row, column) =
          0.03 * std::sin(static_cast<double>((row + 2) * (column + 1)));
    }
  }
  const Eigen::MatrixXd exterior = explicit_exterior(determinants, orbital);
  const xmvb::vb::ExteriorOrbitalTransform transform(determinants, orbital);

  Eigen::MatrixXd coefficients(exterior.rows(), width);
  for (int column = 0; column < coefficients.cols(); ++column) {
    for (int row = 0; row < coefficients.rows(); ++row) {
      coefficients(row, column) =
          std::sin(0.17 * static_cast<double>((row + 1) * (column + 2)));
    }
  }
  Eigen::MatrixXd left = coefficients;
  transform.apply_left(&left);
  require(
      (left - exterior * coefficients).cwiseAbs().maxCoeff() < 2.0e-13,
      "left exterior action differs from explicit compound matrix");

  Eigen::MatrixXd right(width - 2, exterior.rows());
  for (int column = 0; column < right.cols(); ++column) {
    for (int row = 0; row < right.rows(); ++row) {
      right(row, column) =
          std::cos(0.11 * static_cast<double>((row + 3) * (column + 1)));
    }
  }
  const Eigen::MatrixXd original_right = right;
  transform.apply_right(&right);
  require(
      (right - original_right * exterior.transpose())
              .cwiseAbs()
              .maxCoeff() < 2.0e-13,
      "right exterior action differs from explicit compound matrix");

  Eigen::MatrixXd block_storage = Eigen::MatrixXd::Zero(
      original_right.rows(), original_right.cols() + 4);
  block_storage.middleCols(2, original_right.cols()) = original_right;
  auto right_block =
      block_storage.middleCols(2, original_right.cols());
  transform.apply_right_block(right_block);
  require(
      (right_block - right).cwiseAbs().maxCoeff() < 2.0e-13 &&
          block_storage.leftCols(2).isZero(0.0) &&
          block_storage.rightCols(2).isZero(0.0),
      "right exterior block action differs from the owning-matrix path");

  transform.apply_adjoint_right(&right);
  require(
      (right - original_right * exterior.transpose() * exterior)
              .cwiseAbs()
              .maxCoeff() < 4.0e-13,
      "right exterior adjoint is inconsistent");
  transform.apply_adjoint_right_block(right_block);
  require(
      (right_block - right).cwiseAbs().maxCoeff() < 4.0e-13 &&
          block_storage.leftCols(2).isZero(0.0) &&
          block_storage.rightCols(2).isZero(0.0),
      "right exterior block adjoint differs from the owning-matrix path");

  Eigen::MatrixXd adjoint = coefficients;
  transform.apply_adjoint_left(&adjoint);
  require(
      (adjoint - exterior.transpose() * coefficients)
              .cwiseAbs()
              .maxCoeff() < 2.0e-13,
      "left exterior adjoint differs from explicit compound matrix");
  require(transform.shear_pair_count() > 0, "test transform has no shears");

  transform.apply_inverse_left(&left);
  require((left - coefficients).cwiseAbs().maxCoeff() < 3.0e-13,
          "inverse exterior left action fails a forward round trip");
  transform.apply_inverse_adjoint_left(&adjoint);
  require((adjoint - coefficients).cwiseAbs().maxCoeff() < 3.0e-13,
          "inverse exterior adjoint fails a forward round trip");
  Eigen::MatrixXd inverse = coefficients;
  transform.apply_inverse_left(&inverse);
  require((exterior * inverse - coefficients).cwiseAbs().maxCoeff() < 3.0e-13,
          "inverse exterior differs from an independent compound-matrix solve");
  Eigen::MatrixXd inverse_adjoint = coefficients;
  transform.apply_inverse_adjoint_left(&inverse_adjoint);
  require((exterior.transpose() * inverse_adjoint - coefficients).cwiseAbs().maxCoeff() < 3.0e-13,
          "inverse exterior transpose differs from an independent compound solve");
  Eigen::MatrixXd inverse_right = original_right;
  transform.apply_right(&inverse_right);
  transform.apply_inverse_right_block(inverse_right);
  require((inverse_right - original_right).cwiseAbs().maxCoeff() < 3.0e-13,
          "inverse exterior right action fails a forward round trip");
  transform.apply_adjoint_right(&inverse_right);
  transform.apply_inverse_adjoint_right_block(inverse_right);
  require((inverse_right - original_right).cwiseAbs().maxCoeff() < 3.0e-13,
          "inverse exterior right adjoint fails a forward round trip");
  block_storage.middleCols(2, original_right.cols()) = original_right;
  transform.apply_inverse_right_block(right_block);
  const Eigen::MatrixXd right_inverse_reference =
      exterior.partialPivLu().solve(original_right.transpose()).transpose();
  require((right_block - right_inverse_reference).cwiseAbs().maxCoeff() < 3.0e-13 &&
              block_storage.leftCols(2).isZero(0.0) &&
              block_storage.rightCols(2).isZero(0.0),
          "inverse exterior block differs from a compound solve or modifies neighbors");
  block_storage.middleCols(2, original_right.cols()) = original_right;
  transform.apply_inverse_adjoint_right_block(right_block);
  const Eigen::MatrixXd right_inverse_adjoint_reference =
      exterior.transpose().partialPivLu().solve(original_right.transpose()).transpose();
  require((right_block - right_inverse_adjoint_reference).cwiseAbs().maxCoeff() < 3.0e-13 &&
              block_storage.leftCols(2).isZero(0.0) &&
              block_storage.rightCols(2).isZero(0.0),
          "inverse exterior block adjoint differs from a compound solve or modifies neighbors");

  // Non-owning row blocks have an outer stride larger than their row count.
  // Check both their result and untouched neighboring storage with parallel
  // widths that do not divide evenly among workers.
  Eigen::MatrixXd strided = Eigen::MatrixXd::Zero(right.rows() + 4, right.cols());
  strided.middleRows(2, right.rows()) = original_right;
  auto rows = strided.middleRows(2, right.rows());
  transform.apply_right_block(rows);
  require((rows - original_right * exterior.transpose()).cwiseAbs().maxCoeff() < 2.0e-13 &&
              strided.topRows(2).isZero(0.0) && strided.bottomRows(2).isZero(0.0),
          "parallel exterior transform modifies neighbors of a strided block");
  transform.apply_inverse_right_block(rows);
  require((rows - original_right).cwiseAbs().maxCoeff() < 3.0e-13,
          "strided exterior inverse does not recover its input");
}

void check_planner() {
  const auto cerras_alpha = complete_space(11, 5);
  const auto cerras_beta = complete_space(11, 6);
  const auto cerras = xmvb::vb::plan_orthogonal_direct_ci_action(
      cerras_alpha,
      cerras_beta,
      11,
      1);
  require(
      cerras.carrier_complete(),
      "CERRAS-size fixed-spin carrier was not recognized");
  require(cerras.n_alpha_strings == 462 && cerras.n_beta_strings == 462,
          "CERRAS-size determinant count is wrong");
  require(cerras.alpha_same_spin_connections == 181 &&
              cerras.beta_same_spin_connections == 181 &&
              cerras.opposite_spin_connections == 900,
          "CERRAS-size connection count is wrong");
  require(cerras.favors_direct_ci(), "CERRAS-size FLOP plan rejected direct CI");

  const auto loflea_space = complete_space(12, 6);
  const auto loflea = xmvb::vb::plan_orthogonal_direct_ci_action(
      loflea_space,
      loflea_space,
      12,
      1);
  require(
      loflea.carrier_complete(),
      "LOFLEA-size fixed-spin carrier was not recognized");
  require(loflea.n_alpha_strings == 924 && loflea.n_beta_strings == 924,
          "LOFLEA-size determinant count is wrong");
  require(loflea.alpha_same_spin_connections == 262 &&
              loflea.opposite_spin_connections == 1296,
          "LOFLEA-size connection count is wrong");
  require(loflea.favors_direct_ci(), "LOFLEA-size FLOP plan rejected direct CI");

  auto incomplete = cerras_alpha;
  incomplete.pop_back();
  const auto rejected = xmvb::vb::plan_orthogonal_direct_ci_action(
      incomplete,
      cerras_beta,
      11,
      1);
  require(!rejected.carrier_complete() && !rejected.favors_direct_ci(),
          "incomplete fixed-spin space was admitted");

  auto duplicate = cerras_alpha;
  duplicate.back() = duplicate.front();
  const auto duplicate_rejected = xmvb::vb::plan_orthogonal_direct_ci_action(
      duplicate,
      cerras_beta,
      11,
      1);
  require(
      !duplicate_rejected.carrier_complete(),
      "count-complete carrier with a duplicate string was admitted");

  constexpr std::size_t workspace_budget = 256ULL * 1024ULL * 1024ULL;
  require(
      loflea.fits_action_workspace(workspace_budget),
      "LOFLEA-size direct-CI carrier exceeded the action workspace");
  const auto large_half_filled = complete_space(16, 8);
  const auto large_plan = xmvb::vb::plan_orthogonal_direct_ci_action(
      large_half_filled,
      large_half_filled,
      16,
      1);
  require(
      large_plan.carrier_complete(),
      "large fixed-spin carrier was not recognized");
  require(
      !large_plan.fits_action_workspace(workspace_budget),
      "multi-gigabyte direct-CI carrier passed the action workspace limit");
}

void check_sigma_action() {
  constexpr int n_orbitals = 4;
  const auto alpha = complete_space(n_orbitals, 2);
  const auto beta = complete_space(n_orbitals, 2);
  Eigen::MatrixXd orbital = Eigen::MatrixXd::Zero(n_orbitals, n_orbitals);
  orbital <<
      1.10, 0.08, -0.03, 0.05,
      0.00, 0.91,  0.06, 0.02,
      0.00, 0.00,  1.07, 0.04,
      0.00, 0.00,  0.00, 0.96;
  const Eigen::MatrixXd overlap = orbital.transpose() * orbital;
  std::vector<double> packed_overlap(
      overlap.data(), overlap.data() + overlap.size());

  Eigen::MatrixXd one_electron(n_orbitals, n_orbitals);
  for (int column = 0; column < n_orbitals; ++column) {
    for (int row = 0; row <= column; ++row) {
      const double value =
          0.13 * std::cos(static_cast<double>((row + 2) * (column + 1)));
      one_electron(row, column) = value;
      one_electron(column, row) = value;
    }
  }
  const int n_pairs = xmvb::vb::packed_active_pair_count(n_orbitals);
  Eigen::MatrixXd pair_factor(7, n_pairs);
  for (int column = 0; column < n_pairs; ++column) {
    for (int row = 0; row < pair_factor.rows(); ++row) {
      pair_factor(row, column) =
          0.07 * std::sin(static_cast<double>((row + 1) * (column + 3)));
    }
  }
  const Eigen::MatrixXd pair_kernel = pair_factor.transpose() * pair_factor;
  xmvb::vb::ActiveSpaceTwoElectronResult two_electron;
  two_electron.packed_active_two_electron_integrals.resize(
      xmvb::vb::packed_active_two_electron_integral_count(n_orbitals));
  for (int column = 0; column < n_pairs; ++column) {
    for (int row = 0; row <= column; ++row) {
      two_electron.packed_active_two_electron_integrals[
          xmvb::vb::TwoElectronIndexer::packed_pair_of_pairs_index(
              row, column)] = pair_kernel(row, column);
    }
  }

  const auto orthogonal = xmvb::vb::orthogonalize_active_integrals(
      packed_overlap,
      one_electron,
      two_electron,
      n_orbitals);

  Eigen::MatrixXd overlap_direction(n_orbitals, n_orbitals);
  Eigen::MatrixXd one_electron_direction(n_orbitals, n_orbitals);
  for (int column = 0; column < n_orbitals; ++column) {
    for (int row = 0; row <= column; ++row) {
      const double overlap_value =
          0.02 * std::sin(static_cast<double>((row + 1) * (column + 3)));
      const double one_electron_value =
          0.03 * std::cos(static_cast<double>((row + 2) * (column + 4)));
      overlap_direction(row, column) = overlap_value;
      overlap_direction(column, row) = overlap_value;
      one_electron_direction(row, column) = one_electron_value;
      one_electron_direction(column, row) = one_electron_value;
    }
  }
  std::vector<double> packed_overlap_direction(
      overlap_direction.data(),
      overlap_direction.data() + overlap_direction.size());
  std::vector<double> packed_two_electron_direction(
      two_electron.packed_active_two_electron_integrals.size());
  for (std::size_t index = 0;
       index < packed_two_electron_direction.size();
       ++index) {
    packed_two_electron_direction[index] =
        0.01 * std::sin(0.13 * static_cast<double>(index + 1));
  }
  const auto orthogonal_direction =
      xmvb::vb::orthogonalize_active_integral_direction(
          orthogonal,
          packed_overlap_direction,
          one_electron_direction,
          packed_two_electron_direction,
          n_orbitals);
  constexpr double direction_step = 1.0e-6;
  std::vector<double> plus_overlap = packed_overlap;
  std::vector<double> minus_overlap = packed_overlap;
  auto plus_two_electron = two_electron;
  auto minus_two_electron = two_electron;
  for (std::size_t index = 0; index < plus_overlap.size(); ++index) {
    plus_overlap[index] += direction_step * packed_overlap_direction[index];
    minus_overlap[index] -= direction_step * packed_overlap_direction[index];
  }
  for (std::size_t index = 0;
       index < packed_two_electron_direction.size();
       ++index) {
    plus_two_electron.packed_active_two_electron_integrals[index] +=
        direction_step * packed_two_electron_direction[index];
    minus_two_electron.packed_active_two_electron_integrals[index] -=
        direction_step * packed_two_electron_direction[index];
  }
  const auto plus_orthogonal = xmvb::vb::orthogonalize_active_integrals(
      plus_overlap,
      one_electron + direction_step * one_electron_direction,
      plus_two_electron,
      n_orbitals);
  const auto minus_orthogonal = xmvb::vb::orthogonalize_active_integrals(
      minus_overlap,
      one_electron - direction_step * one_electron_direction,
      minus_two_electron,
      n_orbitals);
  require(
      (((plus_orthogonal.orbital_transform -
         minus_orthogonal.orbital_transform) /
        (2.0 * direction_step)) -
       orthogonal_direction.orbital_transform).cwiseAbs().maxCoeff() < 2.0e-9,
      "orthogonal orbital-transform direction fails finite differences");
  require(
      (((plus_orthogonal.one_electron - minus_orthogonal.one_electron) /
        (2.0 * direction_step)) -
       orthogonal_direction.one_electron).cwiseAbs().maxCoeff() < 2.0e-9,
      "orthogonal one-electron direction fails finite differences");
  require(
      (((plus_orthogonal.pair_kernel - minus_orthogonal.pair_kernel) /
        (2.0 * direction_step)) -
       orthogonal_direction.pair_kernel).cwiseAbs().maxCoeff() < 2.0e-9,
      "orthogonal two-electron direction fails finite differences");

  const xmvb::vb::ExteriorOrbitalTransform alpha_transform(
      alpha, orthogonal.orbital_transform);
  const xmvb::vb::ExteriorOrbitalTransform beta_transform(
      beta, orthogonal.orbital_transform);
  const xmvb::vb::DirectCiSigmaAction sigma_action(
      alpha, beta, orthogonal);
  const xmvb::vb::DirectCiSigmaAction sigma_direction(
      alpha, beta, orthogonal_direction);

  const Eigen::MatrixXd diagonal=sigma_action.diagonal(orthogonal.one_electron);
  for (Eigen::Index a=0;a<diagonal.rows();++a) for (Eigen::Index b=0;b<diagonal.cols();++b) {
    Eigen::MatrixXd probe=Eigen::MatrixXd::Zero(diagonal.rows(),diagonal.cols());
    probe(a,b)=1.0;
    require(std::abs(sigma_action.apply(probe)(a,b)-diagonal(a,b))<1e-12,
        "occupation-only carrier diagonal differs from exact sigma diagonal");
  }

  constexpr int block_width = 2;
  Eigen::MatrixXd coefficients(
      static_cast<int>(alpha.size()),
      block_width * static_cast<int>(beta.size()));
  for (int column = 0; column < coefficients.cols(); ++column) {
    for (int row = 0; row < coefficients.rows(); ++row) {
      coefficients(row, column) =
          std::sin(0.19 * static_cast<double>((row + 1) * (column + 2)));
    }
  }

  const xmvb::vb::DirectCiSigmaAction shared_direction(
      sigma_action, orthogonal_direction);
  require((shared_direction.apply(coefficients) - sigma_direction.apply(coefficients))
              .cwiseAbs().maxCoeff() < 2.0e-12,
          "shared spin-string links changed the directional sigma action");
  const auto surviving_direction = [&] {
    const xmvb::vb::DirectCiSigmaAction temporary(alpha, beta, orthogonal);
    return xmvb::vb::DirectCiSigmaAction(temporary, orthogonal_direction);
  }();
  require((surviving_direction.apply(coefficients) - sigma_direction.apply(coefficients))
              .cwiseAbs().maxCoeff() < 2.0e-12,
          "shared spin-string links did not survive their original action");
  const xmvb::vb::DirectCiSigmaAction equal_spin_action(alpha, alpha, orthogonal);
  const xmvb::vb::DirectCiSigmaAction equal_spin_direction(
      equal_spin_action, orthogonal_direction);
  const xmvb::vb::DirectCiSigmaAction equal_spin_reference(
      alpha, alpha, orthogonal_direction);
  const Eigen::MatrixXd equal_spin_values = Eigen::MatrixXd::Ones(
      alpha.size(), 2 * alpha.size());
  require((equal_spin_direction.apply(equal_spin_values) -
           equal_spin_reference.apply(equal_spin_values)).cwiseAbs().maxCoeff() < 2.0e-12,
          "shared alpha/beta links changed the equal-spin sigma action");

  Eigen::MatrixXd left_coefficients(coefficients.rows(), coefficients.cols());
  for (int column = 0; column < left_coefficients.cols(); ++column) {
    for (int row = 0; row < left_coefficients.rows(); ++row) {
      left_coefficients(row, column) =
          std::cos(0.23 * static_cast<double>((row + 2) * (column + 1)));
    }
  }
  const xmvb::vb::DirectCiIntegralAdjoint integral_adjoint =
      sigma_action.integral_adjoint(left_coefficients, coefficients);
  check_fused_response(sigma_action, left_coefficients, coefficients);
  check_fused_response(sigma_action, coefficients, coefficients);
  check_fused_response(sigma_action, Eigen::MatrixXd::Zero(
      coefficients.rows(), coefficients.cols()), coefficients);
  Eigen::MatrixXd padded_left = Eigen::MatrixXd::Zero(
      coefficients.rows() + 3, coefficients.cols());
  Eigen::MatrixXd padded_right = padded_left;
  padded_left.middleRows(1, coefficients.rows()) = left_coefficients;
  padded_right.middleRows(2, coefficients.rows()) = coefficients;
  check_fused_response(sigma_action, padded_left.middleRows(1, coefficients.rows()),
                       padded_right.middleRows(2, coefficients.rows()));
  const xmvb::vb::DirectCiIntegralAdjoint reverse_integral_adjoint =
      sigma_action.integral_adjoint(coefficients, left_coefficients);
  require(
      (integral_adjoint.one_electron -
       reverse_integral_adjoint.one_electron).cwiseAbs().maxCoeff() < 2.0e-12,
      "direct-CI one-electron adjoint violates bilinear reciprocity");
  require(
      (integral_adjoint.pair_kernel -
       reverse_integral_adjoint.pair_kernel).cwiseAbs().maxCoeff() < 2.0e-12,
      "direct-CI pair-kernel adjoint violates bilinear reciprocity");
  const double predicted_integral_direction =
      (integral_adjoint.one_electron.cwiseProduct(
           orthogonal_direction.one_electron)).sum() +
      (integral_adjoint.pair_kernel.cwiseProduct(
           orthogonal_direction.pair_kernel)).sum();
  const double applied_integral_direction =
      (left_coefficients.cwiseProduct(
           sigma_direction.apply(coefficients))).sum();
  require(
      std::abs(predicted_integral_direction - applied_integral_direction) <
          3.0e-10 * std::max(1.0, std::abs(applied_integral_direction)),
      "direct-CI integral adjoint fails the bilinear directional identity");
  const auto pack_pair_kernel = [n_orbitals](
      const Eigen::MatrixXd& kernel,
      xmvb::vb::ActiveSpaceTwoElectronResult* result) {
    result->packed_active_two_electron_integrals.resize(
        xmvb::vb::packed_active_two_electron_integral_count(n_orbitals));
    for (int column = 0; column < kernel.cols(); ++column) {
      for (int row = 0; row <= column; ++row) {
        result->packed_active_two_electron_integrals[
            xmvb::vb::TwoElectronIndexer::packed_pair_of_pairs_index(
                row, column)] = kernel(row, column);
      }
    }
  };
  constexpr double adjoint_step = 1.0e-6;
  auto plus_sigma_integrals = orthogonal;
  auto minus_sigma_integrals = orthogonal;
  plus_sigma_integrals.one_electron +=
      adjoint_step * orthogonal_direction.one_electron;
  minus_sigma_integrals.one_electron -=
      adjoint_step * orthogonal_direction.one_electron;
  plus_sigma_integrals.pair_kernel +=
      adjoint_step * orthogonal_direction.pair_kernel;
  minus_sigma_integrals.pair_kernel -=
      adjoint_step * orthogonal_direction.pair_kernel;
  pack_pair_kernel(
      plus_sigma_integrals.pair_kernel,
      &plus_sigma_integrals.two_electron);
  pack_pair_kernel(
      minus_sigma_integrals.pair_kernel,
      &minus_sigma_integrals.two_electron);
  const double plus_bilinear =
      (left_coefficients.cwiseProduct(
           xmvb::vb::DirectCiSigmaAction(
               alpha, beta, plus_sigma_integrals).apply(coefficients))).sum();
  const double minus_bilinear =
      (left_coefficients.cwiseProduct(
           xmvb::vb::DirectCiSigmaAction(
               alpha, beta, minus_sigma_integrals).apply(coefficients))).sum();
  const double finite_difference_bilinear =
      (plus_bilinear - minus_bilinear) / (2.0 * adjoint_step);
  require(
      std::abs(predicted_integral_direction - finite_difference_bilinear) <
          2.0e-8 * std::max(1.0, std::abs(finite_difference_bilinear)),
      "direct-CI integral adjoint fails central finite differences");

  const auto alpha_transform_direction = alpha_transform.direction(
      orthogonal_direction.orbital_transform);
  Eigen::MatrixXd transformed_coefficients = coefficients;
  Eigen::MatrixXd transformed_direction = Eigen::MatrixXd::Zero(
      coefficients.rows(), coefficients.cols());
  alpha_transform.apply_directional_left(
      alpha_transform_direction,
      &transformed_coefficients,
      &transformed_direction);
  Eigen::MatrixXd plus_transformed = coefficients;
  Eigen::MatrixXd minus_transformed = coefficients;
  xmvb::vb::ExteriorOrbitalTransform(
      alpha, plus_orthogonal.orbital_transform).apply_left(&plus_transformed);
  xmvb::vb::ExteriorOrbitalTransform(
      alpha, minus_orthogonal.orbital_transform).apply_left(&minus_transformed);
  require(
      (((plus_transformed - minus_transformed) /
        (2.0 * direction_step)) -
       transformed_direction).cwiseAbs().maxCoeff() < 2.0e-9,
      "exterior-transform direction fails finite differences");

  Eigen::MatrixXd right_probe(5, static_cast<int>(alpha.size()));
  for (int column = 0; column < right_probe.cols(); ++column) {
    for (int row = 0; row < right_probe.rows(); ++row) {
      right_probe(row, column) =
          std::cos(0.09 * static_cast<double>((row + 2) * (column + 1)));
    }
  }
  Eigen::MatrixXd directional_right = right_probe;
  Eigen::MatrixXd delta_directional_right = Eigen::MatrixXd::Zero(
      right_probe.rows(), right_probe.cols());
  alpha_transform.apply_directional_right(
      alpha_transform_direction,
      &directional_right,
      &delta_directional_right);
  Eigen::MatrixXd plus_right = right_probe;
  Eigen::MatrixXd minus_right = right_probe;
  xmvb::vb::ExteriorOrbitalTransform(
      alpha, plus_orthogonal.orbital_transform).apply_right(&plus_right);
  xmvb::vb::ExteriorOrbitalTransform(
      alpha, minus_orthogonal.orbital_transform).apply_right(&minus_right);
  require(
      (((plus_right - minus_right) / (2.0 * direction_step)) -
       delta_directional_right).cwiseAbs().maxCoeff() < 2.0e-9,
      "right exterior-transform direction fails finite differences");

  Eigen::MatrixXd directional_adjoint_left = coefficients;
  Eigen::MatrixXd delta_directional_adjoint_left = Eigen::MatrixXd::Zero(
      coefficients.rows(), coefficients.cols());
  alpha_transform.apply_directional_adjoint_left(
      alpha_transform_direction,
      &directional_adjoint_left,
      &delta_directional_adjoint_left);
  Eigen::MatrixXd plus_adjoint_left = coefficients;
  Eigen::MatrixXd minus_adjoint_left = coefficients;
  xmvb::vb::ExteriorOrbitalTransform(
      alpha,
      plus_orthogonal.orbital_transform).apply_adjoint_left(
          &plus_adjoint_left);
  xmvb::vb::ExteriorOrbitalTransform(
      alpha,
      minus_orthogonal.orbital_transform).apply_adjoint_left(
          &minus_adjoint_left);
  require(
      (((plus_adjoint_left - minus_adjoint_left) /
        (2.0 * direction_step)) -
       delta_directional_adjoint_left).cwiseAbs().maxCoeff() < 2.0e-9,
      "adjoint-left exterior direction fails finite differences");

  Eigen::MatrixXd directional_adjoint_right = right_probe;
  Eigen::MatrixXd delta_directional_adjoint_right = Eigen::MatrixXd::Zero(
      right_probe.rows(), right_probe.cols());
  alpha_transform.apply_directional_adjoint_right(
      alpha_transform_direction,
      &directional_adjoint_right,
      &delta_directional_adjoint_right);
  Eigen::MatrixXd plus_adjoint_right = right_probe;
  Eigen::MatrixXd minus_adjoint_right = right_probe;
  xmvb::vb::ExteriorOrbitalTransform(
      alpha,
      plus_orthogonal.orbital_transform).apply_adjoint_right(
          &plus_adjoint_right);
  xmvb::vb::ExteriorOrbitalTransform(
      alpha,
      minus_orthogonal.orbital_transform).apply_adjoint_right(
          &minus_adjoint_right);
  require(
      (((plus_adjoint_right - minus_adjoint_right) /
        (2.0 * direction_step)) -
       delta_directional_adjoint_right).cwiseAbs().maxCoeff() < 2.0e-9,
      "adjoint-right exterior direction fails finite differences");

  Eigen::MatrixXd sigma_probe = coefficients;
  const Eigen::MatrixXd plus_sigma =
      xmvb::vb::DirectCiSigmaAction(
          alpha, beta, plus_orthogonal).apply(sigma_probe);
  const Eigen::MatrixXd minus_sigma =
      xmvb::vb::DirectCiSigmaAction(
          alpha, beta, minus_orthogonal).apply(sigma_probe);
  require(
      (((plus_sigma - minus_sigma) / (2.0 * direction_step)) -
       sigma_direction.apply(sigma_probe)).cwiseAbs().maxCoeff() < 2.0e-9,
      "direct-CI sigma direction fails finite differences");

  const auto apply_orthogonal_action = [&alpha, &beta, block_width](
      const xmvb::vb::OrthogonalActiveIntegrals& integrals,
      const Eigen::MatrixXd& input) {
    const xmvb::vb::ExteriorOrbitalTransform alpha_action_transform(
        alpha, integrals.orbital_transform);
    const xmvb::vb::ExteriorOrbitalTransform beta_action_transform(
        beta, integrals.orbital_transform);
    const xmvb::vb::DirectCiSigmaAction action_sigma(
        alpha, beta, integrals);
    Eigen::MatrixXd values = input;
    alpha_action_transform.apply_left(&values);
    for (int block = 0; block < block_width; ++block) {
      Eigen::MatrixXd value_block = values.middleCols(
          block * static_cast<int>(beta.size()),
          static_cast<int>(beta.size()));
      beta_action_transform.apply_right(&value_block);
      values.middleCols(
          block * static_cast<int>(beta.size()),
          static_cast<int>(beta.size())) = value_block;
    }
    Eigen::MatrixXd hamiltonian = action_sigma.apply(values);
    Eigen::MatrixXd overlap_image = values;
    for (int block = 0; block < block_width; ++block) {
      Eigen::MatrixXd hamiltonian_block = hamiltonian.middleCols(
          block * static_cast<int>(beta.size()),
          static_cast<int>(beta.size()));
      Eigen::MatrixXd overlap_block = overlap_image.middleCols(
          block * static_cast<int>(beta.size()),
          static_cast<int>(beta.size()));
      beta_action_transform.apply_adjoint_right(&hamiltonian_block);
      beta_action_transform.apply_adjoint_right(&overlap_block);
      hamiltonian.middleCols(
          block * static_cast<int>(beta.size()),
          static_cast<int>(beta.size())) = hamiltonian_block;
      overlap_image.middleCols(
          block * static_cast<int>(beta.size()),
          static_cast<int>(beta.size())) = overlap_block;
    }
    alpha_action_transform.apply_adjoint_left(&hamiltonian);
    alpha_action_transform.apply_adjoint_left(&overlap_image);
    return std::make_pair(hamiltonian, overlap_image);
  };

  Eigen::MatrixXd action_values = coefficients;
  Eigen::MatrixXd delta_action_values = Eigen::MatrixXd::Zero(
      coefficients.rows(), coefficients.cols());
  alpha_transform.apply_directional_left(
      alpha_transform_direction,
      &action_values,
      &delta_action_values);
  const auto beta_transform_direction = beta_transform.direction(
      orthogonal_direction.orbital_transform);
  for (int block = 0; block < block_width; ++block) {
    Eigen::MatrixXd value_block = action_values.middleCols(
        block * static_cast<int>(beta.size()),
        static_cast<int>(beta.size()));
    Eigen::MatrixXd delta_value_block = delta_action_values.middleCols(
        block * static_cast<int>(beta.size()),
        static_cast<int>(beta.size()));
    beta_transform.apply_directional_right(
        beta_transform_direction,
        &value_block,
        &delta_value_block);
    action_values.middleCols(
        block * static_cast<int>(beta.size()),
        static_cast<int>(beta.size())) = value_block;
    delta_action_values.middleCols(
        block * static_cast<int>(beta.size()),
        static_cast<int>(beta.size())) = delta_value_block;
  }
  const Eigen::MatrixXd inverse_orbital_transform =
      orthogonal.orbital_transform
          .template triangularView<Eigen::Upper>()
          .solve(Eigen::MatrixXd::Identity(n_orbitals, n_orbitals));
  const Eigen::MatrixXd relative_orbital_direction =
      orthogonal_direction.orbital_transform * inverse_orbital_transform;
  const Eigen::MatrixXd generator_adjoint =
      sigma_action.one_body_generator_adjoint(
          left_coefficients,
          action_values);
  const double predicted_exterior_direction =
      (generator_adjoint.cwiseProduct(relative_orbital_direction)).sum();
  const double applied_exterior_direction =
      (left_coefficients.cwiseProduct(delta_action_values)).sum();
  require(
      std::abs(predicted_exterior_direction - applied_exterior_direction) <
          3.0e-10 * std::max(1.0, std::abs(applied_exterior_direction)),
      "direct-CI orbital-generator adjoint is inconsistent: predicted=" +
          std::to_string(predicted_exterior_direction) +
          ", applied=" + std::to_string(applied_exterior_direction));
  const Eigen::MatrixXd generator_action =
      sigma_action.apply_one_body_generator(
          action_values,
          relative_orbital_direction);
  require(
      (generator_action - delta_action_values).cwiseAbs().maxCoeff() < 2.0e-11,
      "direct-CI one-body generator action is inconsistent with the "
      "exterior derivative");
  Eigen::MatrixXd batched_forward, batched_transpose;
  const Eigen::MatrixXd generator_transpose = relative_orbital_direction.transpose();
  sigma_action.apply_one_body_generators({
      {action_values, relative_orbital_direction, &batched_forward},
      {left_coefficients, generator_transpose, &batched_transpose}});
  require((batched_forward - generator_action).cwiseAbs().maxCoeff() < 2.0e-12,
          "batched generator changed the exterior derivative");
  const Eigen::MatrixXd expected_transpose =
      sigma_action.apply_one_body_generator(left_coefficients, generator_transpose);
  require((batched_transpose - expected_transpose).cwiseAbs().maxCoeff() < 2.0e-12,
          "batched generator mixed coefficient or generator views");
  Eigen::MatrixXd action_hamiltonian = sigma_action.apply(action_values);
  Eigen::MatrixXd delta_action_hamiltonian =
      sigma_action.apply(delta_action_values) +
      sigma_direction.apply(action_values);
  Eigen::MatrixXd action_overlap = action_values;
  Eigen::MatrixXd delta_action_overlap = delta_action_values;
  for (int block = 0; block < block_width; ++block) {
    Eigen::MatrixXd hamiltonian_block = action_hamiltonian.middleCols(
        block * static_cast<int>(beta.size()),
        static_cast<int>(beta.size()));
    Eigen::MatrixXd delta_hamiltonian_block =
        delta_action_hamiltonian.middleCols(
            block * static_cast<int>(beta.size()),
            static_cast<int>(beta.size()));
    Eigen::MatrixXd overlap_block = action_overlap.middleCols(
        block * static_cast<int>(beta.size()),
        static_cast<int>(beta.size()));
    Eigen::MatrixXd delta_overlap_block =
        delta_action_overlap.middleCols(
            block * static_cast<int>(beta.size()),
            static_cast<int>(beta.size()));
    beta_transform.apply_directional_adjoint_right(
        beta_transform_direction,
        &hamiltonian_block,
        &delta_hamiltonian_block);
    beta_transform.apply_directional_adjoint_right(
        beta_transform_direction,
        &overlap_block,
        &delta_overlap_block);
    action_hamiltonian.middleCols(
        block * static_cast<int>(beta.size()),
        static_cast<int>(beta.size())) = hamiltonian_block;
    delta_action_hamiltonian.middleCols(
        block * static_cast<int>(beta.size()),
        static_cast<int>(beta.size())) = delta_hamiltonian_block;
    action_overlap.middleCols(
        block * static_cast<int>(beta.size()),
        static_cast<int>(beta.size())) = overlap_block;
    delta_action_overlap.middleCols(
        block * static_cast<int>(beta.size()),
        static_cast<int>(beta.size())) = delta_overlap_block;
  }
  alpha_transform.apply_directional_adjoint_left(
      alpha_transform_direction,
      &action_hamiltonian,
      &delta_action_hamiltonian);
  alpha_transform.apply_directional_adjoint_left(
      alpha_transform_direction,
      &action_overlap,
      &delta_action_overlap);
  const auto plus_action = apply_orthogonal_action(
      plus_orthogonal, coefficients);
  const auto minus_action = apply_orthogonal_action(
      minus_orthogonal, coefficients);
  require(
      ((((plus_action.first - minus_action.first) /
         (2.0 * direction_step)) -
        delta_action_hamiltonian).cwiseAbs().maxCoeff()) < 5.0e-9,
      "complete direct-CI Hamiltonian direction fails finite differences");
  require(
      ((((plus_action.second - minus_action.second) /
         (2.0 * direction_step)) -
        delta_action_overlap).cwiseAbs().maxCoeff()) < 5.0e-9,
      "complete direct-CI overlap direction fails finite differences");

  constexpr double test_energy = 0.37;
  const Eigen::MatrixXd state_coefficients = coefficients.leftCols(
      static_cast<int>(beta.size()));
  const Eigen::MatrixXd orthogonal_state = action_values.leftCols(
      static_cast<int>(beta.size()));
  const Eigen::MatrixXd orthogonal_sigma = sigma_action.apply(
      orthogonal_state);
  const Eigen::MatrixXd orthogonal_residual =
      orthogonal_sigma - test_energy * orthogonal_state;
  const xmvb::vb::DirectCiIntegralAdjoint state_integral_adjoint =
      sigma_action.integral_adjoint(
          orthogonal_state,
          orthogonal_state);
  const Eigen::MatrixXd state_generator_adjoint =
      2.0 * sigma_action.one_body_generator_adjoint(
          orthogonal_residual,
          orthogonal_state);
  const xmvb::vb::NonorthogonalActiveIntegralAdjoint state_adjoint =
      xmvb::vb::backpropagate_orthogonal_active_integral_adjoint(
          orthogonal,
          state_integral_adjoint.one_electron,
          state_integral_adjoint.pair_kernel,
          state_generator_adjoint);
  Eigen::MatrixXd original_pair_direction(n_pairs, n_pairs);
  for (int column = 0; column < n_pairs; ++column) {
    for (int row = 0; row < n_pairs; ++row) {
      original_pair_direction(row, column) =
          packed_two_electron_direction[
              xmvb::vb::TwoElectronIndexer::packed_pair_of_pairs_index(
                  row, column)];
    }
  }
  const double predicted_state_direction =
      (state_adjoint.overlap.cwiseProduct(overlap_direction)).sum() +
      (state_adjoint.one_electron.cwiseProduct(
           one_electron_direction)).sum() +
      (state_adjoint.pair_kernel.cwiseProduct(
           original_pair_direction)).sum();
  const auto state_lagrangian = [test_energy, &state_coefficients](
      const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& images) {
    return (state_coefficients.cwiseProduct(
        images.first.leftCols(state_coefficients.cols()) -
        test_energy * images.second.leftCols(
            state_coefficients.cols()))).sum();
  };
  const double finite_difference_state_direction =
      (state_lagrangian(plus_action) - state_lagrangian(minus_action)) /
      (2.0 * direction_step);
  require(
      std::abs(predicted_state_direction -
               finite_difference_state_direction) <
          3.0e-8 *
          std::max(1.0, std::abs(finite_difference_state_direction)),
      "direct-CI active-integral adjoint fails finite differences");

  const Eigen::MatrixXd response_coefficients =
      left_coefficients.leftCols(static_cast<int>(beta.size()));
  Eigen::MatrixXd orthogonal_response = response_coefficients;
  alpha_transform.apply_left(&orthogonal_response);
  beta_transform.apply_right(&orthogonal_response);
  const Eigen::MatrixXd response_residual =
      sigma_action.apply(orthogonal_response) -
      test_energy * orthogonal_response;
  const xmvb::vb::DirectCiIntegralAdjoint response_integral_adjoint =
      sigma_action.integral_adjoint(
          orthogonal_response, orthogonal_state);
  const Eigen::MatrixXd response_generator_adjoint = 2.0 *
      (sigma_action.one_body_generator_adjoint(
           response_residual, orthogonal_state) +
       sigma_action.one_body_generator_adjoint(
           orthogonal_residual, orthogonal_response));
  const xmvb::vb::NonorthogonalActiveIntegralAdjoint response_adjoint =
      xmvb::vb::backpropagate_orthogonal_active_integral_adjoint(
          orthogonal,
          2.0 * response_integral_adjoint.one_electron,
          2.0 * response_integral_adjoint.pair_kernel,
          response_generator_adjoint);
  const double predicted_response_direction =
      (response_adjoint.overlap.cwiseProduct(overlap_direction)).sum() +
      (response_adjoint.one_electron.cwiseProduct(
           one_electron_direction)).sum() +
      (response_adjoint.pair_kernel.cwiseProduct(
           original_pair_direction)).sum();
  const double applied_response_direction = 2.0 *
      (response_coefficients.cwiseProduct(
           delta_action_hamiltonian.leftCols(response_coefficients.cols()) -
           test_energy *
               delta_action_overlap.leftCols(response_coefficients.cols())))
          .sum();
  require(
      std::abs(predicted_response_direction - applied_response_direction) <
          3.0e-10 * std::max(1.0, std::abs(applied_response_direction)),
      "direct-CI response adjoint is inconsistent: predicted=" +
          std::to_string(predicted_response_direction) +
          ", applied=" + std::to_string(applied_response_direction));

  const Eigen::MatrixXd directional_orthogonal_state =
      delta_action_values.leftCols(static_cast<int>(beta.size()));
  const Eigen::MatrixXd directional_orthogonal_sigma =
      sigma_action.apply(directional_orthogonal_state) +
      sigma_direction.apply(orthogonal_state);
  const Eigen::MatrixXd directional_orthogonal_residual =
      directional_orthogonal_sigma -
      test_energy * directional_orthogonal_state;
  const xmvb::vb::DirectCiIntegralAdjoint left_state_adjoint_direction =
      sigma_action.integral_adjoint(
          directional_orthogonal_state,
          orthogonal_state);
  const xmvb::vb::DirectCiIntegralAdjoint right_state_adjoint_direction =
      sigma_action.integral_adjoint(
          orthogonal_state,
          directional_orthogonal_state);
  const Eigen::MatrixXd one_gradient_direction =
      left_state_adjoint_direction.one_electron +
      right_state_adjoint_direction.one_electron;
  const Eigen::MatrixXd pair_gradient_direction =
      left_state_adjoint_direction.pair_kernel +
      right_state_adjoint_direction.pair_kernel;
  const Eigen::MatrixXd generator_gradient_direction = 2.0 *
      (sigma_action.one_body_generator_adjoint(
           directional_orthogonal_residual,
           orthogonal_state) +
       sigma_action.one_body_generator_adjoint(
           orthogonal_residual,
           directional_orthogonal_state));
  const xmvb::vb::NonorthogonalActiveIntegralAdjoint
      analytic_state_adjoint_direction =
          xmvb::vb::backpropagate_orthogonal_active_integral_adjoint_direction(
              orthogonal,
              orthogonal_direction,
              state_integral_adjoint.one_electron,
              state_integral_adjoint.pair_kernel,
              state_generator_adjoint,
              one_gradient_direction,
              pair_gradient_direction,
              generator_gradient_direction);
  const auto build_state_adjoint = [
      &alpha,
      &beta,
      &state_coefficients,
      test_energy](const xmvb::vb::OrthogonalActiveIntegrals& integrals) {
    const xmvb::vb::ExteriorOrbitalTransform alpha_state_transform(
        alpha, integrals.orbital_transform);
    const xmvb::vb::ExteriorOrbitalTransform beta_state_transform(
        beta, integrals.orbital_transform);
    const xmvb::vb::DirectCiSigmaAction state_sigma(alpha, beta, integrals);
    Eigen::MatrixXd state = state_coefficients;
    alpha_state_transform.apply_left(&state);
    beta_state_transform.apply_right(&state);
    const Eigen::MatrixXd sigma = state_sigma.apply(state);
    const Eigen::MatrixXd residual = sigma - test_energy * state;
    const xmvb::vb::DirectCiIntegralAdjoint integral_gradient =
        state_sigma.integral_adjoint(state, state);
    const Eigen::MatrixXd generator_gradient =
        2.0 * state_sigma.one_body_generator_adjoint(residual, state);
    return xmvb::vb::backpropagate_orthogonal_active_integral_adjoint(
        integrals,
        integral_gradient.one_electron,
        integral_gradient.pair_kernel,
        generator_gradient);
  };
  const xmvb::vb::NonorthogonalActiveIntegralAdjoint plus_state_adjoint =
      build_state_adjoint(plus_orthogonal);
  const xmvb::vb::NonorthogonalActiveIntegralAdjoint minus_state_adjoint =
      build_state_adjoint(minus_orthogonal);
  require(
      ((((plus_state_adjoint.overlap - minus_state_adjoint.overlap) /
         (2.0 * direction_step)) -
        analytic_state_adjoint_direction.overlap).cwiseAbs().maxCoeff()) <
          2.0e-7,
      "direct-CI overlap-adjoint direction fails finite differences");
  require(
      ((((plus_state_adjoint.one_electron -
          minus_state_adjoint.one_electron) /
         (2.0 * direction_step)) -
        analytic_state_adjoint_direction.one_electron)
           .cwiseAbs()
           .maxCoeff()) < 2.0e-7,
      "direct-CI one-electron-adjoint direction fails finite differences");
  require(
      ((((plus_state_adjoint.pair_kernel -
          minus_state_adjoint.pair_kernel) /
         (2.0 * direction_step)) -
        analytic_state_adjoint_direction.pair_kernel)
           .cwiseAbs()
           .maxCoeff()) < 2.0e-7,
      "direct-CI pair-kernel-adjoint direction fails finite differences");

  Eigen::MatrixXd orthogonal_coefficients = coefficients;
  for (int block = 0; block < block_width; ++block) {
    Eigen::MatrixXd block_coefficients =
        orthogonal_coefficients.middleCols(
            block * static_cast<int>(beta.size()),
            static_cast<int>(beta.size()));
    alpha_transform.apply_left(&block_coefficients);
    beta_transform.apply_right(&block_coefficients);
    orthogonal_coefficients.middleCols(
        block * static_cast<int>(beta.size()),
        static_cast<int>(beta.size())) = block_coefficients;
  }
  Eigen::MatrixXd direct_hamiltonian =
      sigma_action.apply(orthogonal_coefficients);
  Eigen::MatrixXd direct_overlap = orthogonal_coefficients;
  for (int block = 0; block < block_width; ++block) {
    Eigen::MatrixXd hamiltonian_block = direct_hamiltonian.middleCols(
        block * static_cast<int>(beta.size()),
        static_cast<int>(beta.size()));
    Eigen::MatrixXd overlap_block = direct_overlap.middleCols(
        block * static_cast<int>(beta.size()),
        static_cast<int>(beta.size()));
    beta_transform.apply_adjoint_right(&hamiltonian_block);
    alpha_transform.apply_adjoint_left(&hamiltonian_block);
    beta_transform.apply_adjoint_right(&overlap_block);
    alpha_transform.apply_adjoint_left(&overlap_block);
    direct_hamiltonian.middleCols(
        block * static_cast<int>(beta.size()),
        static_cast<int>(beta.size())) = hamiltonian_block;
    direct_overlap.middleCols(
        block * static_cast<int>(beta.size()),
        static_cast<int>(beta.size())) = overlap_block;
  }

  const int product_dimension =
      static_cast<int>(alpha.size() * beta.size());
  Eigen::MatrixXd dense_hamiltonian(product_dimension, product_dimension);
  Eigen::MatrixXd dense_overlap(product_dimension, product_dimension);
  const xmvb::vb::DeterminantPairEvaluator evaluator;
  for (int target_alpha = 0;
       target_alpha < static_cast<int>(alpha.size());
       ++target_alpha) {
    for (int target_beta = 0;
         target_beta < static_cast<int>(beta.size());
         ++target_beta) {
      const int target =
          target_alpha * static_cast<int>(beta.size()) + target_beta;
      for (int source_alpha = 0;
           source_alpha < static_cast<int>(alpha.size());
           ++source_alpha) {
        for (int source_beta = 0;
             source_beta < static_cast<int>(beta.size());
             ++source_beta) {
          const int source =
              source_alpha * static_cast<int>(beta.size()) + source_beta;
          const auto pair = evaluator.evaluate(
              alpha[target_alpha],
              alpha[source_alpha],
              beta[target_beta],
              beta[source_beta],
              packed_overlap,
              one_electron,
              n_orbitals,
              two_electron,
              false);
          dense_hamiltonian(target, source) = pair.total_hamiltonian;
          dense_overlap(target, source) = pair.overlap_determinant;
        }
      }
    }
  }

  Eigen::MatrixXd coefficient_columns(product_dimension, block_width);
  for (int block = 0; block < block_width; ++block) {
    for (int alpha_index = 0;
         alpha_index < static_cast<int>(alpha.size());
         ++alpha_index) {
      for (int beta_index = 0;
           beta_index < static_cast<int>(beta.size());
           ++beta_index) {
        coefficient_columns(
            alpha_index * static_cast<int>(beta.size()) + beta_index,
            block) = coefficients(
                alpha_index,
                block * static_cast<int>(beta.size()) + beta_index);
      }
    }
  }
  const Eigen::MatrixXd reference_hamiltonian =
      dense_hamiltonian * coefficient_columns;
  const Eigen::MatrixXd reference_overlap =
      dense_overlap * coefficient_columns;
  for (int block = 0; block < block_width; ++block) {
    for (int alpha_index = 0;
         alpha_index < static_cast<int>(alpha.size());
         ++alpha_index) {
      for (int beta_index = 0;
           beta_index < static_cast<int>(beta.size());
           ++beta_index) {
        const int product =
            alpha_index * static_cast<int>(beta.size()) + beta_index;
        const int column =
            block * static_cast<int>(beta.size()) + beta_index;
        require(
            std::abs(
                direct_hamiltonian(alpha_index, column) -
                reference_hamiltonian(product, block)) < 2.0e-10,
            "orthogonal direct-CI Hamiltonian action is inconsistent");
        require(
            std::abs(
                direct_overlap(alpha_index, column) -
                reference_overlap(product, block)) < 2.0e-10,
            "orthogonal direct-CI overlap action is inconsistent");
      }
    }
  }

  // Project the exact determinant-product action onto three arbitrary sparse
  // structure vectors. The structure subspace is deliberately tiny and is
  // not invariant under H; exactness therefore cannot rely on a full-CI
  // structure basis.
  constexpr int n_test_structures = 3;
  Eigen::MatrixXd structure_expansion = Eigen::MatrixXd::Zero(
      product_dimension, n_test_structures);
  for (int structure = 0; structure < n_test_structures; ++structure) {
    for (int term = 0; term < 7; ++term) {
      const int product =
          (5 * term + 11 * structure + structure * term) % product_dimension;
      structure_expansion(product, structure) +=
          std::sin(0.31 * static_cast<double>((term + 1) * (structure + 2)));
    }
  }
  Eigen::MatrixXd structure_vectors(n_test_structures, block_width);
  for (int column = 0; column < block_width; ++column) {
    for (int row = 0; row < n_test_structures; ++row) {
      structure_vectors(row, column) =
          std::cos(0.27 * static_cast<double>((row + 1) * (column + 3)));
    }
  }
  const Eigen::MatrixXd product_vectors =
      structure_expansion * structure_vectors;
  Eigen::MatrixXd packed_product_vectors(
      static_cast<int>(alpha.size()),
      block_width * static_cast<int>(beta.size()));
  for (int block = 0; block < block_width; ++block) {
    for (int alpha_index = 0;
         alpha_index < static_cast<int>(alpha.size());
         ++alpha_index) {
      for (int beta_index = 0;
           beta_index < static_cast<int>(beta.size());
           ++beta_index) {
        packed_product_vectors(
            alpha_index,
            block * static_cast<int>(beta.size()) + beta_index) =
            product_vectors(
                alpha_index * static_cast<int>(beta.size()) + beta_index,
                block);
      }
    }
  }
  const auto projected_action = apply_orthogonal_action(
      orthogonal, packed_product_vectors);
  Eigen::MatrixXd projected_hamiltonian_products(
      product_dimension, block_width);
  Eigen::MatrixXd projected_overlap_products(
      product_dimension, block_width);
  for (int block = 0; block < block_width; ++block) {
    for (int alpha_index = 0;
         alpha_index < static_cast<int>(alpha.size());
         ++alpha_index) {
      for (int beta_index = 0;
           beta_index < static_cast<int>(beta.size());
           ++beta_index) {
        const int product =
            alpha_index * static_cast<int>(beta.size()) + beta_index;
        const int packed_column =
            block * static_cast<int>(beta.size()) + beta_index;
        projected_hamiltonian_products(product, block) =
            projected_action.first(alpha_index, packed_column);
        projected_overlap_products(product, block) =
            projected_action.second(alpha_index, packed_column);
      }
    }
  }
  const Eigen::MatrixXd direct_three_structure_hamiltonian =
      structure_expansion.transpose() * projected_hamiltonian_products;
  const Eigen::MatrixXd direct_three_structure_overlap =
      structure_expansion.transpose() * projected_overlap_products;
  const Eigen::MatrixXd reference_three_structure_hamiltonian =
      structure_expansion.transpose() * dense_hamiltonian *
      structure_expansion * structure_vectors;
  const Eigen::MatrixXd reference_three_structure_overlap =
      structure_expansion.transpose() * dense_overlap *
      structure_expansion * structure_vectors;
  require(
      (direct_three_structure_hamiltonian -
       reference_three_structure_hamiltonian).cwiseAbs().maxCoeff() < 2.0e-10,
      "direct CI is inexact on an arbitrary three-structure subspace");
  require(
      (direct_three_structure_overlap -
       reference_three_structure_overlap).cwiseAbs().maxCoeff() < 2.0e-10,
      "direct-CI overlap is inexact on an arbitrary three-structure subspace");

  std::vector<std::vector<int>> product_alpha;
  std::vector<std::vector<int>> product_beta;
  std::vector<std::vector<xmvb::vb::StructureExpansionTerm>>
      determinant_to_structure_terms(product_dimension);
  product_alpha.reserve(product_dimension);
  product_beta.reserve(product_dimension);
  for (int alpha_index = 0;
       alpha_index < static_cast<int>(alpha.size());
       ++alpha_index) {
    for (int beta_index = 0;
         beta_index < static_cast<int>(beta.size());
         ++beta_index) {
      const int product =
          alpha_index * static_cast<int>(beta.size()) + beta_index;
      product_alpha.push_back(alpha[alpha_index]);
      product_beta.push_back(beta[beta_index]);
      for (int structure = 0; structure < n_test_structures; ++structure) {
        const double coefficient = structure_expansion(product, structure);
        if (coefficient != 0.0) {
          determinant_to_structure_terms[product].push_back(
              xmvb::vb::StructureExpansionTerm{structure, coefficient});
        }
      }
    }
  }
  const auto topology = xmvb::vb::build_same_spin_pair_topology(
      product_alpha, product_beta, n_orbitals);
  require(!topology.enabled(), "topology-only context materialized pair data");
  require(topology.alpha_pair_cache.empty() && topology.beta_pair_cache.empty(),
          "topology-only context retained ordered pair arrays");

  xmvb::vb::StructureDiagonal exact_diagonal;
  exact_diagonal.hamiltonian =
      (structure_expansion.transpose() * dense_hamiltonian *
       structure_expansion).diagonal();
  exact_diagonal.overlap =
      (structure_expansion.transpose() * dense_overlap *
       structure_expansion).diagonal();
  const xmvb::vb::FullDeterminantStructureHamiltonianOverlapBuilder builder;
  const xmvb::vb::StructureDiagonal streamed_diagonal =
      builder.build_diagonal(
          product_alpha,
          product_beta,
          determinant_to_structure_terms,
          packed_overlap,
          one_electron,
          n_orbitals,
          two_electron,
          n_test_structures,
          topology);
  require(
      (streamed_diagonal.hamiltonian - exact_diagonal.hamiltonian)
              .cwiseAbs()
              .maxCoeff() < 2.0e-10,
      "streamed structure Hamiltonian diagonal is inexact");
  require(
      (streamed_diagonal.overlap - exact_diagonal.overlap)
              .cwiseAbs()
              .maxCoeff() < 2.0e-10,
      "streamed structure overlap diagonal is inexact");
  const Eigen::VectorXd overlap_only_diagonal =
      builder.build_exact_overlap_diagonal(
          product_alpha,
          product_beta,
          determinant_to_structure_terms,
          packed_overlap,
          n_orbitals,
          n_test_structures,
          topology);
  require(
      (overlap_only_diagonal - exact_diagonal.overlap)
              .cwiseAbs()
              .maxCoeff() < 2.0e-10,
      "overlap-only structure diagonal is inexact");

  xmvb::vb::StructureDiagonal expected_jacobi;
  expected_jacobi.hamiltonian =
      Eigen::VectorXd::Zero(n_test_structures);
  expected_jacobi.overlap = Eigen::VectorXd::Zero(n_test_structures);
  for (int determinant = 0; determinant < product_dimension; ++determinant) {
    for (int structure = 0; structure < n_test_structures; ++structure) {
      const double coefficient =
          structure_expansion(determinant, structure);
      expected_jacobi.hamiltonian[structure] +=
          coefficient * coefficient *
          dense_hamiltonian(determinant, determinant);
      expected_jacobi.overlap[structure] +=
          coefficient * coefficient * dense_overlap(determinant, determinant);
    }
  }
  const xmvb::vb::StructureDiagonal jacobi =
      builder.build_davidson_jacobi_preconditioner(
          product_alpha,
          product_beta,
          determinant_to_structure_terms,
          packed_overlap,
          one_electron,
          n_orbitals,
          two_electron,
          n_test_structures,
          topology);
  require(
      (jacobi.hamiltonian - expected_jacobi.hamiltonian)
              .cwiseAbs()
              .maxCoeff() < 2.0e-10,
      "determinant-diagonal Hamiltonian Jacobi model is inexact");
  require(
      (jacobi.overlap - expected_jacobi.overlap)
              .cwiseAbs()
              .maxCoeff() < 2.0e-10,
      "determinant-diagonal overlap Jacobi model is inexact");
  require(
      (jacobi.hamiltonian - exact_diagonal.hamiltonian)
              .cwiseAbs()
              .maxCoeff() > 1.0e-8,
      "Jacobi model unexpectedly retained determinant cross terms");
  xmvb::vb::StructureAction topology_action(
      determinant_to_structure_terms,
      n_test_structures,
      topology,
      packed_overlap,
      one_electron,
      two_electron,
      n_orbitals,
      &exact_diagonal);
  require(
      topology_action.storage().orthogonal_direct_ci,
      "complete fixed-spin carrier did not select direct CI");
  bool selected_metric_rejected = false;
  try {
    const xmvb::vb::StructureMetricInverse invalid(
        topology_action, topology.alpha_reuse_table.unique_determinants,
        topology.beta_reuse_table.unique_determinants, false, 65536, 65536);
  } catch (const std::invalid_argument&) {
    selected_metric_rejected = true;
  }
  require(selected_metric_rejected,
          "complete carrier incorrectly admitted a selected structure metric inverse");

  // Complete spin-product basis with nonorthogonal recoupling inside each
  // spatial occupation pattern. This tests nontrivial Gram blocks and chunks.
  std::vector<std::vector<xmvb::vb::StructureExpansionTerm>> full_terms(product_dimension);
  Eigen::MatrixXd full_expansion = Eigen::MatrixXd::Identity(product_dimension, product_dimension);
  std::map<std::vector<int>, int> first_product;
  for (int product = 0; product < product_dimension; ++product) {
    full_terms[product].push_back({product, 1.0});
    std::vector<int> occupation(n_orbitals, 0);
    for (int orbital : product_alpha[product]) ++occupation[orbital];
    for (int orbital : product_beta[product]) ++occupation[orbital];
    const auto [entry, inserted] = first_product.emplace(occupation, product);
    if (!inserted) {
      full_terms[entry->second].push_back({product, 0.25});
      full_expansion(entry->second, product) = 0.25;
    }
  }
  const Eigen::MatrixXd full_overlap = full_expansion.transpose() * dense_overlap * full_expansion;
  xmvb::vb::StructureDiagonal full_diagonal;
  full_diagonal.overlap = full_overlap.diagonal();
  full_diagonal.hamiltonian =
      (full_expansion.transpose() * dense_hamiltonian * full_expansion).diagonal();
  const xmvb::vb::StructureAction full_action(
      full_terms, product_dimension, topology, packed_overlap,
      one_electron, two_electron, n_orbitals, &full_diagonal);
  check_joint_structure_response(full_action, alpha.size(), beta.size());
  check_joint_structure_response(topology_action, alpha.size(), beta.size());
  const std::size_t one_column_bytes = 5 * product_dimension * sizeof(double);
  require(xmvb::vb::StructureMetricInverse::fits(
      full_action, topology.alpha_reuse_table.unique_determinants,
      topology.beta_reuse_table.unique_determinants, 65536, one_column_bytes),
      "structure metric planner rejected an affordable block factorization");
  const xmvb::vb::StructureMetricInverse full_metric(
      full_action, topology.alpha_reuse_table.unique_determinants,
      topology.beta_reuse_table.unique_determinants, true, 65536, one_column_bytes);
  Eigen::MatrixXd metric_rhs(product_dimension, 7);
  for (int column = 0; column < metric_rhs.cols(); ++column)
    for (int row = 0; row < metric_rhs.rows(); ++row)
      metric_rhs(row, column) = std::sin(0.19 * (row + 1) * (column + 2));
  const Eigen::MatrixXd metric_solution = full_metric.apply(metric_rhs);
  for (int column = 0; column < metric_rhs.cols(); ++column) {
    const Eigen::VectorXd coefficients = metric_rhs.col(column);
    const double energy = 0.2 * (column + 1);
    const auto images = full_action.apply(coefficients);
    const Eigen::MatrixXd residual = images.hamiltonian - energy * images.overlap;
    const double reference = std::sqrt(
        residual.col(0).dot(full_metric.apply(residual).col(0)) /
        coefficients.dot(images.overlap.col(0)));
    require(std::abs(full_action.complete_space_residual_norm(coefficients, energy) - reference) <
                1e-10 * std::max(1.0, reference),
            "orthogonal carrier residual differs from the physical dual metric");
  }
  const Eigen::MatrixXd reference_metric_solution = full_overlap.partialPivLu().solve(metric_rhs);
  require((metric_solution - reference_metric_solution).cwiseAbs().maxCoeff() < 3.0e-11 &&
              (full_action.apply(metric_solution).overlap - metric_rhs).cwiseAbs().maxCoeff() < 3.0e-11,
          "sparse complete-space metric inverse differs from its dense reference");
  require(full_metric.factor_bytes() < static_cast<std::size_t>(product_dimension) * product_dimension * sizeof(double) &&
              full_metric.max_block_size() == 6 && full_metric.block_count() == 19,
          "spin metric lost its occupation-block representation");
  // Independent explicit pullback of the carrier diagonal, including BOTH
  // Riesz maps. A single F^T D^-1 F is not a structure-coordinate inverse.
  const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(product_dimension,product_dimension);
  const Eigen::MatrixXd carriers = full_action.orthogonalize_structure_block(identity);
  Eigen::MatrixXd f(product_dimension,product_dimension);
  for (int column=0;column<product_dimension;++column)
    for (int a=0;a<static_cast<int>(alpha.size());++a)
      for (int b=0;b<static_cast<int>(beta.size());++b)
        f(a*beta.size()+b,column)=carriers(a,column*beta.size()+b);
  const Eigen::MatrixXd inverse_m = full_overlap.partialPivLu().solve(identity);
  const double root_energy = diagonal.minCoeff()-1.0;
  Eigen::VectorXd inverse_d(product_dimension);
  for (int a=0;a<static_cast<int>(alpha.size());++a)
    for (int b=0;b<static_cast<int>(beta.size());++b)
      inverse_d[a*beta.size()+b]=1.0/(diagonal(a,b)-root_energy);
  const Eigen::MatrixXd reference_p =
      inverse_m*f.transpose()*inverse_d.asDiagonal()*f*inverse_m;
  const Eigen::MatrixXd riesz = full_metric.apply(metric_rhs.leftCols(1));
  const Eigen::MatrixXd carrier_image = full_action.apply_inverse_carrier_diagonal(
      riesz,Eigen::VectorXd::Constant(1,root_energy),0.0);
  const Eigen::MatrixXd correction = full_metric.apply(carrier_image);
  const Eigen::MatrixXd expected_correction = reference_p*metric_rhs.leftCols(1);
  require((correction-expected_correction).norm()<
              1e-10*std::max(1.0,expected_correction.norm()),
          "carrier Davidson inverse differs from the dense metric pullback");
  auto original_action = std::make_unique<xmvb::vb::StructureAction>(
      full_terms,product_dimension,topology,packed_overlap,
      one_electron,two_electron,n_orbitals,&full_diagonal);
  xmvb::vb::StructureMetricInverse moved_metric(
      *original_action,topology.alpha_reuse_table.unique_determinants,
      topology.beta_reuse_table.unique_determinants,true,65536,one_column_bytes);
  xmvb::vb::StructureAction moved_action(std::move(*original_action));
  moved_metric.rebind(moved_action);
  original_action.reset();
  require((moved_metric.apply(metric_rhs)-metric_solution).norm()<1e-11,
          "moving a forward action invalidated its reusable metric factors");
  for (bool factor_limit : {true, false}) {
    require(!xmvb::vb::StructureMetricInverse::fits(
        full_action, topology.alpha_reuse_table.unique_determinants,
        topology.beta_reuse_table.unique_determinants,
        factor_limit ? 0 : 65536,
        factor_limit ? one_column_bytes : one_column_bytes - 1),
        "structure metric planner ignored its byte budget");
    bool rejected = false;
    try {
      const xmvb::vb::StructureMetricInverse invalid(
          full_action, topology.alpha_reuse_table.unique_determinants,
          topology.beta_reuse_table.unique_determinants, true,
          factor_limit ? 0 : 65536,
          factor_limit ? one_column_bytes : one_column_bytes - 1);
    } catch (const std::length_error&) {
      rejected = true;
    }
    require(rejected, "structure metric ignored its factor or workspace byte limit");
  }
  const auto topology_images = topology_action.apply(structure_vectors);
  const xmvb::vb::StructureAction action_only(
      determinant_to_structure_terms, n_test_structures, topology,
      packed_overlap, one_electron, two_electron, n_orbitals, nullptr,
      xmvb::vb::StructureActionMode::ActionOnly);
  const auto bare_images = action_only.apply(structure_vectors);
  require(action_only.storage().diagonal_bytes == 0 &&
              action_only.preconditioner_diagonal().overlap.size() == 0 &&
              (bare_images.hamiltonian - topology_images.hamiltonian).norm() < 2e-10 &&
              (bare_images.overlap - topology_images.overlap).norm() < 2e-10,
          "action-only direct CI built a diagonal or changed H/S");
  Eigen::MatrixXd normalized = structure_vectors;
  auto scaled = bare_images;
  for (int col = 0; col < normalized.cols(); ++col) {
    const double norm = std::sqrt(normalized.col(col).dot(scaled.overlap.col(col)));
    normalized.col(col) /= norm;
    scaled.hamiltonian.col(col) /= norm;
    scaled.overlap.col(col) /= norm;
  }
  const auto normalized_images = action_only.apply(normalized);
  require((scaled.hamiltonian - normalized_images.hamiltonian).norm() < 2e-10 &&
              (scaled.overlap - normalized_images.overlap).norm() < 2e-10,
          "column normalization requires a second H/S action");
  require((topology_action.apply_overlap(structure_vectors) -
      topology_images.overlap).cwiseAbs().maxCoeff() < 2.0e-12,
      "overlap-only direct-CI action differs from the combined action");
  require(
      (topology_images.hamiltonian -
       reference_three_structure_hamiltonian).cwiseAbs().maxCoeff() < 2.0e-10,
      "topology-only Hamiltonian action differs from the dense reference");
  require(
      (topology_images.overlap -
       reference_three_structure_overlap).cwiseAbs().maxCoeff() < 2.0e-10,
      "topology-only overlap action differs from the dense reference");
  std::vector<std::vector<int>> incomplete_alpha;
  std::vector<std::vector<int>> incomplete_beta;
  std::vector<std::vector<xmvb::vb::StructureExpansionTerm>>
      incomplete_terms;
  std::vector<int> retained_products;
  for (int product = 0; product < product_dimension; ++product) {
    if (product_alpha[product] == alpha.back()) {
      continue;
    }
    incomplete_alpha.push_back(product_alpha[product]);
    incomplete_beta.push_back(product_beta[product]);
    incomplete_terms.push_back(determinant_to_structure_terms[product]);
    retained_products.push_back(product);
  }
  Eigen::MatrixXd incomplete_expansion(
      static_cast<int>(retained_products.size()), n_test_structures);
  Eigen::MatrixXd incomplete_hamiltonian(
      retained_products.size(), retained_products.size());
  Eigen::MatrixXd incomplete_overlap(
      retained_products.size(), retained_products.size());
  for (int left = 0; left < static_cast<int>(retained_products.size()); ++left) {
    incomplete_expansion.row(left) =
        structure_expansion.row(retained_products[left]);
    for (int right = 0;
         right < static_cast<int>(retained_products.size());
         ++right) {
      incomplete_hamiltonian(left, right) =
          dense_hamiltonian(retained_products[left], retained_products[right]);
      incomplete_overlap(left, right) =
          dense_overlap(retained_products[left], retained_products[right]);
    }
  }
  const auto incomplete_topology = xmvb::vb::build_same_spin_pair_topology(
      incomplete_alpha, incomplete_beta, n_orbitals);
  xmvb::vb::StructureAction streamed_action(
      incomplete_terms,
      n_test_structures,
      incomplete_topology,
      packed_overlap,
      one_electron,
      two_electron,
      n_orbitals,
      &exact_diagonal);
  require(
      !streamed_action.storage().orthogonal_direct_ci,
      "incomplete fixed-spin carrier incorrectly selected direct CI");
  const auto streamed_images =
      streamed_action.apply_streamed(structure_vectors);
  const xmvb::vb::StructureAction pair_action_only(
      incomplete_terms, n_test_structures, incomplete_topology,
      packed_overlap, one_electron, two_electron, n_orbitals, nullptr,
      xmvb::vb::StructureActionMode::ActionOnly);
  const auto bare_pair_images = pair_action_only.apply(structure_vectors);
  require(pair_action_only.storage().diagonal_bytes == 0 &&
              (bare_pair_images.hamiltonian - streamed_images.hamiltonian).norm() < 2e-10 &&
              (bare_pair_images.overlap - streamed_images.overlap).norm() < 2e-10,
          "action-only incomplete carrier built a diagonal or changed H/S");
  require((streamed_action.apply_overlap(structure_vectors) -
      streamed_images.overlap).cwiseAbs().maxCoeff() < 2.0e-12,
      "overlap-only pair action differs from the streamed action");
  const Eigen::MatrixXd incomplete_reference_hamiltonian =
      incomplete_expansion.transpose() * incomplete_hamiltonian *
      incomplete_expansion * structure_vectors;
  const Eigen::MatrixXd incomplete_reference_overlap =
      incomplete_expansion.transpose() * incomplete_overlap *
      incomplete_expansion * structure_vectors;
  require(
      (streamed_images.hamiltonian -
       incomplete_reference_hamiltonian).cwiseAbs().maxCoeff() < 2.0e-10,
      "streamed Hamiltonian action differs from the dense reference");
  require(
      (streamed_images.overlap -
       incomplete_reference_overlap).cwiseAbs().maxCoeff() < 2.0e-10,
      "streamed overlap action differs from the dense reference");
}

/** @brief Exercise partial final tiles and multiple RHS with a known spectrum. */
void check_sigma_tile_boundary() {
  constexpr int n_orbitals = 10;
  constexpr int block_width = 2;
  // 252 strings: partial alpha tile and uneven beta ownership at 32 workers.
  const auto spin = complete_space(n_orbitals, 5);
  xmvb::vb::OrthogonalActiveIntegrals integrals;
  integrals.one_electron = Eigen::MatrixXd::Zero(n_orbitals, n_orbitals);
  for (int i = 0; i < n_orbitals; ++i)
    integrals.one_electron(i, i) = 0.17 * (i + 1);
  const int n_pairs = xmvb::vb::packed_active_pair_count(n_orbitals);
  integrals.pair_kernel = Eigen::MatrixXd::Zero(n_pairs, n_pairs);
  const xmvb::vb::DirectCiSigmaAction action(spin, spin, integrals);
  const int n_strings = static_cast<int>(spin.size());
  Eigen::MatrixXd coefficients(n_strings, block_width * n_strings);
  Eigen::VectorXd energies = Eigen::VectorXd::Zero(n_strings);
  for (int i = 0; i < n_strings; ++i)
    for (int orbital : spin[i]) energies[i] += integrals.one_electron(orbital, orbital);
  Eigen::MatrixXd expected(coefficients.rows(), coefficients.cols());
  for (int j = 0; j < coefficients.cols(); ++j)
    for (int i = 0; i < coefficients.rows(); ++i) {
      coefficients(i, j) = std::sin(0.031 * (i + 1) * (j + 2));
      expected(i, j) = (energies[i] + energies[j % n_strings]) * coefficients(i, j);
    }
  require((action.apply(coefficients) - expected).cwiseAbs().maxCoeff() < 2e-12,
          "direct-CI partial tile or RHS packing differs from known spectrum");
  check_fused_response(action, coefficients, coefficients);
#ifdef _OPENMP
  const int dynamic = omp_get_dynamic();
  omp_set_dynamic(1);
  const Eigen::MatrixXd dynamic_result = action.apply(coefficients);
  check_fused_response(action, coefficients, coefficients);
  omp_set_dynamic(dynamic);
  require((dynamic_result - expected).cwiseAbs().maxCoeff() < 2e-12,
          "direct-CI dynamic team reduction differs from known spectrum");
  Eigen::MatrixXd nested_result;
#pragma omp parallel num_threads(2)
  {
#pragma omp single
    {
      nested_result = action.apply(coefficients);
      check_fused_response(action, coefficients, coefficients);
    }
  }
  require((nested_result - expected).cwiseAbs().maxCoeff() < 2e-12,
          "direct-CI nested team reduction differs from known spectrum");
#endif
  const auto adjoint = action.integral_adjoint(coefficients, coefficients);
  require(std::abs(adjoint.one_electron.cwiseProduct(integrals.one_electron).sum() -
                   coefficients.cwiseProduct(expected).sum()) < 2e-8,
          "tiled direct-CI adjoint differs from action bilinear form");
#ifdef _OPENMP
  omp_set_dynamic(1);
  const auto dynamic_adjoint = action.integral_adjoint(coefficients, coefficients);
  omp_set_dynamic(dynamic);
  require((dynamic_adjoint.pair_kernel - adjoint.pair_kernel).norm() <
              1e-11 * std::max(1.0, adjoint.pair_kernel.norm()),
          "adjoint output ownership depends on the actual team size");
  xmvb::vb::DirectCiIntegralAdjoint nested_adjoint;
#pragma omp parallel num_threads(2)
  {
#pragma omp single
    nested_adjoint = action.integral_adjoint(coefficients, coefficients);
  }
  require((nested_adjoint.one_electron - adjoint.one_electron).norm() <
              1e-11 * std::max(1.0, adjoint.one_electron.norm()),
          "nested adjoint differs from the shared-tile contraction");
#endif
}

void check_distinct_spin_sigma_action() {
  constexpr int n_orbitals = 4;
  constexpr int block_width = 2;
  const auto alpha = complete_space(n_orbitals, 1);
  const auto beta = complete_space(n_orbitals, 2);

  Eigen::MatrixXd overlap = Eigen::MatrixXd::Identity(
      n_orbitals, n_orbitals);
  std::vector<double> packed_overlap(
      overlap.data(), overlap.data() + overlap.size());
  Eigen::MatrixXd one_electron(n_orbitals, n_orbitals);
  for (int column = 0; column < n_orbitals; ++column) {
    for (int row = 0; row <= column; ++row) {
      const double value =
          0.11 * std::cos(0.3 * static_cast<double>((row + 1) * (column + 2)));
      one_electron(row, column) = value;
      one_electron(column, row) = value;
    }
  }
  const int n_pairs = xmvb::vb::packed_active_pair_count(n_orbitals);
  Eigen::MatrixXd factor(6, n_pairs);
  for (int column = 0; column < n_pairs; ++column) {
    for (int row = 0; row < factor.rows(); ++row) {
      factor(row, column) =
          0.09 * std::sin(0.17 * static_cast<double>((row + 2) * (column + 1)));
    }
  }
  const Eigen::MatrixXd pair_kernel = factor.transpose() * factor;
  const auto pack_kernel = [n_orbitals](
      const Eigen::MatrixXd& kernel) {
    xmvb::vb::ActiveSpaceTwoElectronResult result;
    result.packed_active_two_electron_integrals.resize(
        xmvb::vb::packed_active_two_electron_integral_count(n_orbitals));
    for (int column = 0; column < kernel.cols(); ++column) {
      for (int row = 0; row <= column; ++row) {
        result.packed_active_two_electron_integrals[
            xmvb::vb::TwoElectronIndexer::packed_pair_of_pairs_index(
                row, column)] = kernel(row, column);
      }
    }
    return result;
  };
  const xmvb::vb::ActiveSpaceTwoElectronResult two_electron =
      pack_kernel(pair_kernel);
  const xmvb::vb::OrthogonalActiveIntegrals integrals{
      Eigen::MatrixXd::Identity(n_orbitals, n_orbitals),
      one_electron,
      pair_kernel,
      two_electron};
  const xmvb::vb::DirectCiSigmaAction action(alpha, beta, integrals);

  Eigen::MatrixXd coefficients(
      static_cast<int>(alpha.size()),
      block_width * static_cast<int>(beta.size()));
  for (int column = 0; column < coefficients.cols(); ++column) {
    for (int row = 0; row < coefficients.rows(); ++row) {
      coefficients(row, column) =
          std::sin(0.21 * static_cast<double>((row + 2) * (column + 1)));
    }
  }
  const int n_beta = static_cast<int>(beta.size());
  const int product_dimension = static_cast<int>(alpha.size()) * n_beta;
  std::vector<std::vector<int>> product_alpha, product_beta;
  std::vector<std::vector<xmvb::vb::StructureExpansionTerm>> terms(product_dimension);
  for (const auto& a : alpha)
    for (const auto& b : beta) {
      const int i = product_alpha.size();
      product_alpha.push_back(a);
      product_beta.push_back(b);
      terms[i].push_back({i, 1.0});
    }
  const auto topology = xmvb::vb::build_same_spin_pair_topology(
      product_alpha, product_beta, n_orbitals);
  const xmvb::vb::StructureAction structure_action(terms, product_dimension,
      topology, packed_overlap, one_electron, two_electron, n_orbitals,
      nullptr, xmvb::vb::StructureActionMode::ActionOnly);
  check_joint_structure_response(structure_action, alpha.size(), beta.size());
  Eigen::MatrixXd dense_hamiltonian(product_dimension, product_dimension);
  const xmvb::vb::DeterminantPairEvaluator evaluator;
  for (int target_alpha = 0;
       target_alpha < static_cast<int>(alpha.size());
       ++target_alpha) {
    for (int target_beta = 0; target_beta < n_beta; ++target_beta) {
      const int target = target_alpha * n_beta + target_beta;
      for (int source_alpha = 0;
           source_alpha < static_cast<int>(alpha.size());
           ++source_alpha) {
        for (int source_beta = 0; source_beta < n_beta; ++source_beta) {
          const int source = source_alpha * n_beta + source_beta;
          dense_hamiltonian(target, source) = evaluator.evaluate(
              alpha[target_alpha],
              alpha[source_alpha],
              beta[target_beta],
              beta[source_beta],
              packed_overlap,
              one_electron,
              n_orbitals,
              two_electron,
              false).total_hamiltonian;
        }
      }
    }
  }
  const Eigen::MatrixXd sigma = action.apply(coefficients);
  check_fused_response(action, coefficients, coefficients);
  for (int block = 0; block < block_width; ++block) {
    Eigen::VectorXd packed_coefficients(product_dimension);
    for (int alpha_index = 0;
         alpha_index < static_cast<int>(alpha.size());
         ++alpha_index) {
      for (int beta_index = 0; beta_index < n_beta; ++beta_index) {
        packed_coefficients[alpha_index * n_beta + beta_index] =
            coefficients(alpha_index, block * n_beta + beta_index);
      }
    }
    const Eigen::VectorXd expected = dense_hamiltonian * packed_coefficients;
    for (int alpha_index = 0;
         alpha_index < static_cast<int>(alpha.size());
         ++alpha_index) {
      for (int beta_index = 0; beta_index < n_beta; ++beta_index) {
        require(
            std::abs(
                sigma(alpha_index, block * n_beta + beta_index) -
                expected[alpha_index * n_beta + beta_index]) < 2.0e-12,
            "distinct-spin direct-CI sigma differs from dense assembly");
      }
    }
  }
}

void check_demand_driven_structure_diagonal() {
  // Put one demanded pair in the first cell of a unique-spin space larger
  // than the forward matrix tile. A tile-backed diagonal lookup would evaluate
  // the entire 256x256 tile even though all other determinants have zero
  // structure coefficient.
  constexpr int n_orbitals = 11;
  constexpr int n_unique_alpha = 257;
  auto alpha = complete_space(n_orbitals, 5);
  alpha.resize(n_unique_alpha);
  std::vector<std::vector<int>> beta(n_unique_alpha, std::vector<int>{0});
  std::vector<std::vector<xmvb::vb::StructureExpansionTerm>> expansion(
      n_unique_alpha);
  expansion.front().push_back({0, 1.0});

  const Eigen::MatrixXd overlap =
      Eigen::MatrixXd::Identity(n_orbitals, n_orbitals);
  const std::vector<double> packed_overlap(
      overlap.data(),
      overlap.data() + overlap.size());
  Eigen::MatrixXd one_electron(n_orbitals, n_orbitals);
  for (int column = 0; column < n_orbitals; ++column) {
    for (int row = 0; row <= column; ++row) {
      const double value =
          0.03 * std::cos(static_cast<double>((row + 1) * (column + 2)));
      one_electron(row, column) = value;
      one_electron(column, row) = value;
    }
  }
  xmvb::vb::ActiveSpaceTwoElectronResult two_electron;
  two_electron.packed_active_two_electron_integrals.assign(
      xmvb::vb::packed_active_two_electron_integral_count(n_orbitals),
      0.0);

  const auto large_topology = xmvb::vb::build_same_spin_pair_topology(
      alpha,
      beta,
      n_orbitals);
  require(!large_topology.enabled(),
          "demand-driven diagonal regression unexpectedly built pair caches");
  require(large_topology.alpha_reuse_table.unique_determinants.size() ==
              n_unique_alpha,
          "demand-driven diagonal regression did not span multiple tiles");

  const xmvb::vb::FullDeterminantStructureHamiltonianOverlapBuilder builder;
  const xmvb::vb::StructureDiagonal large_diagonal = builder.build_diagonal(
      alpha,
      beta,
      expansion,
      packed_overlap,
      one_electron,
      n_orbitals,
      two_electron,
      1,
      large_topology);

  const std::vector<std::vector<int>> reference_alpha{alpha.front()};
  const std::vector<std::vector<int>> reference_beta{beta.front()};
  const std::vector<std::vector<xmvb::vb::StructureExpansionTerm>>
      reference_expansion{{{0, 1.0}}};
  const auto reference_topology = xmvb::vb::build_same_spin_pair_topology(
      reference_alpha,
      reference_beta,
      n_orbitals);
  const xmvb::vb::StructureDiagonal reference_diagonal =
      builder.build_diagonal(
          reference_alpha,
          reference_beta,
          reference_expansion,
          packed_overlap,
          one_electron,
          n_orbitals,
          two_electron,
          1,
          reference_topology);

  require(
      std::abs(
          large_diagonal.hamiltonian[0] -
          reference_diagonal.hamiltonian[0]) < 2.0e-12,
      "unused unique determinants changed the demand-driven Hamiltonian diagonal");
  require(
      std::abs(
          large_diagonal.overlap[0] -
          reference_diagonal.overlap[0]) < 2.0e-12,
      "unused unique determinants changed the demand-driven overlap diagonal");

  std::vector<std::vector<xmvb::vb::StructureExpansionTerm>>
      jacobi_expansion(n_unique_alpha);
  double expected_jacobi_hamiltonian = 0.0;
  double expected_jacobi_overlap = 0.0;
  const xmvb::vb::DeterminantPairEvaluator evaluator;
  for (int determinant = 0; determinant < n_unique_alpha; ++determinant) {
    const double coefficient =
        std::sin(0.013 * static_cast<double>(determinant + 1));
    if (determinant == 0) {
      // Duplicate determinant/structure terms must be combined before the
      // coefficient is squared.
      jacobi_expansion[determinant].push_back({0, 0.4 * coefficient});
      jacobi_expansion[determinant].push_back({0, 0.6 * coefficient});
    } else {
      jacobi_expansion[determinant].push_back({0, coefficient});
    }
    const auto self_pair = evaluator.evaluate(
        alpha[determinant],
        alpha[determinant],
        beta[determinant],
        beta[determinant],
        packed_overlap,
        one_electron,
        n_orbitals,
        two_electron,
        false);
    expected_jacobi_hamiltonian +=
        coefficient * coefficient * self_pair.total_hamiltonian;
    expected_jacobi_overlap +=
        coefficient * coefficient * self_pair.overlap_determinant;
  }
  const xmvb::vb::StructureDiagonal jacobi =
      builder.build_davidson_jacobi_preconditioner(
          alpha,
          beta,
          jacobi_expansion,
          packed_overlap,
          one_electron,
          n_orbitals,
          two_electron,
          1,
          large_topology);
  require(
      std::abs(jacobi.hamiltonian[0] - expected_jacobi_hamiltonian) < 2.0e-10,
      "large determinant-diagonal Hamiltonian Jacobi model is inexact");
  require(
      std::abs(jacobi.overlap[0] - expected_jacobi_overlap) < 2.0e-10,
      "large determinant-diagonal overlap Jacobi model is inexact");

  const std::vector<std::vector<xmvb::vb::StructureExpansionTerm>>
      cancelling_expansion{{{0, 1.0}, {0, -1.0}}};
  bool rejected_zero_overlap = false;
  try {
    static_cast<void>(builder.build_davidson_jacobi_preconditioner(
        reference_alpha,
        reference_beta,
        cancelling_expansion,
        packed_overlap,
        one_electron,
        n_orbitals,
        two_electron,
        1,
        reference_topology));
  } catch (const std::domain_error&) {
    rejected_zero_overlap = true;
  }
  require(
      rejected_zero_overlap,
      "zero Davidson Jacobi overlap was not rejected explicitly");
}

/** @brief Distinct spin sectors, including equal sizes but different occupations. */
void check_fused_spin_sectors() {
  for (const auto& sector : {std::array<int, 3>{11, 5, 6}, {10, 4, 5}, {10, 5, 4},
                            {11, 1, 5}, {11, 5, 1}, {10, 5, 0}, {10, 0, 5}}) {
    const int n_orbitals = sector[0];
    const auto alpha = complete_space(n_orbitals, sector[1]);
    const auto beta = complete_space(n_orbitals, sector[2]);
    const int n_pairs = xmvb::vb::packed_active_pair_count(n_orbitals);
    xmvb::vb::OrthogonalActiveIntegrals integrals;
    integrals.one_electron = Eigen::MatrixXd::Identity(n_orbitals, n_orbitals);
    Eigen::MatrixXd factors(5, n_pairs);
    for (int pair = 0; pair < n_pairs; ++pair)
      for (int row = 0; row < factors.rows(); ++row)
        factors(row, pair) = 0.04 * std::cos(0.13 * (row + 1) * (pair + 2));
    integrals.pair_kernel = factors.transpose() * factors;
    const xmvb::vb::DirectCiSigmaAction action(alpha, beta, integrals);
    Eigen::MatrixXd response(alpha.size(), 2 * beta.size());
    Eigen::MatrixXd state(response.rows(), response.cols());
    for (int column = 0; column < response.cols(); ++column)
      for (int row = 0; row < response.rows(); ++row) {
        response(row, column) = std::sin(0.07 * (row + 2) * (column + 3));
        state(row, column) = std::cos(0.03 * (row + 3) * (column + 1));
      }
    check_fused_response(action, response, state);
  }
}

/** @brief Same-binary joint/separate response comparison, not an optimizer replay. */
void benchmark_structure_response() {
  using Clock = std::chrono::steady_clock;
  std::cout << "kernel,n_orbitals,n_alpha,n_beta,calls,wall_seconds,checksum\n"
            << std::setprecision(17);
  for (const auto sector : {std::array<int,3>{6,3,3}, {11,5,5}, {11,5,6}, {12,6,6}}) {
    const int n = sector[0];
    const auto alpha = complete_space(n, sector[1]), beta = complete_space(n, sector[2]);
    const int products = alpha.size()*beta.size();
    std::vector<std::vector<int>> product_alpha, product_beta;
    std::vector<std::vector<xmvb::vb::StructureExpansionTerm>> terms(products);
    product_alpha.reserve(products);
    product_beta.reserve(products);
    for (const auto& a : alpha)
      for (const auto& b : beta) {
        const int k = product_alpha.size();
        product_alpha.push_back(a);
        product_beta.push_back(b);
        for (int j = 0; j < 3; ++j)
          terms[k].push_back({j, std::sin(0.17*(k+1)*(j+2))/std::sqrt(products)});
      }
    const auto topology = xmvb::vb::build_same_spin_pair_topology(product_alpha, product_beta, n);
    Eigen::MatrixXd r = Eigen::MatrixXd::Identity(n,n);
    for (int i = 0; i < n; ++i)
      for (int j = i+1; j < n; ++j) r(i,j) = 0.02*std::sin((i+1)*(j+2));
    const Eigen::MatrixXd overlap = r.transpose()*r;
    const std::vector<double> packed_overlap(overlap.data(), overlap.data()+overlap.size());
    const Eigen::MatrixXd h = -Eigen::MatrixXd::Identity(n,n);
    const int pairs = xmvb::vb::packed_active_pair_count(n);
    xmvb::vb::ActiveSpaceTwoElectronResult eri;
    eri.packed_active_two_electron_integrals.resize(
        xmvb::vb::packed_active_two_electron_integral_count(n));
    for (int j = 0; j < pairs; ++j)
      for (int i = 0; i <= j; ++i)
        eri.packed_active_two_electron_integrals[
            xmvb::vb::TwoElectronIndexer::packed_pair_of_pairs_index(i,j)] =
            0.03*std::cos(0.07*(i+1)*(j+2));
    const xmvb::vb::StructureAction action(terms, 3, topology, packed_overlap,
        h, eri, n, nullptr, xmvb::vb::StructureActionMode::ActionOnly);
    xmvb::vb::SelectedStateDeterminantMatrices states;
    states.n_structures = 3;
    states.n_unique_alpha = alpha.size();
    states.n_unique_beta = beta.size();
    states.normalized_state_weights = {1.0};
    states.states.resize(1);
    states.states[0].normalized_state_weight = 1.0;
    const Eigen::Vector3d c(0.3,-0.2,0.7), q(0.4,0.6,-0.3);
    states.states[0].coefficient_matrix = action.expand_structure_block(c);
    const auto accepted = action.prepare_active_adjoint(states, {-1.3});
    const auto evaluate = [&](bool joint) {
      xmvb::vb::StructureActionResult images;
      xmvb::vb::StructureResponseSweep sweep;
      Eigen::MatrixXd metric;
      if (joint) {
        sweep = action.apply_response(accepted,q);
        images = std::move(sweep.images);
        metric = images.overlap;
      } else {
        metric = action.apply_overlap(q);
        images = action.apply(q);
      }
      Eigen::MatrixXd multiplier(1,1);
      multiplier(0,0) = -c.dot(images.hamiltonian.col(0)+1.3*images.overlap.col(0));
      xmvb::vb::StructureActiveIntegralAdjoint adjoint;
      if (joint) adjoint = action.finish_response(accepted,std::move(sweep),multiplier);
      else {
        auto response = states;
        response.states[0].coefficient_matrix = action.expand_structure_block(q);
        adjoint = action.active_integral_response_adjoint(accepted,response,multiplier);
      }
      return metric.squaredNorm()+images.hamiltonian.squaredNorm()+
          adjoint.overlap.squaredNorm()+adjoint.one_electron.squaredNorm()+
          adjoint.pair_kernel.squaredNorm();
    };
    const double separate = evaluate(false), joint = evaluate(true);
    require(std::abs(separate-joint) < 3e-11*std::max(1.0,std::abs(separate)),
        "joint response timing checksum differs");
    for (bool fused : {false,true}) {
      constexpr int repeats = 3;
      double checksum = 0;
      const auto start = Clock::now();
      for (int i = 0; i < repeats; ++i) checksum += evaluate(fused);
      std::cout << (fused ? "joint," : "separate,") << n << ',' << alpha.size()
                << ',' << beta.size() << ',' << repeats << ','
                << std::chrono::duration<double>(Clock::now()-start).count()
                << ',' << checksum << '\n';
    }
  }
}

/** @brief Bounded kernel timings; keep separate from convergence benchmarks. */
void benchmark_sigma() {
  using Clock = std::chrono::steady_clock;
  std::cout << "kernel,n_orbitals,n_alpha,n_beta,block_width,calls,wall_seconds,checksum\n"
            << std::setprecision(17);
  for (const auto& sector : {std::array<int, 3>{6, 3, 3}, {11, 5, 5}, {11, 5, 6},
                            {12, 6, 6}, {11, 1, 5}, {11, 5, 1}}) {
    const int n_orbitals = sector[0];
    const auto alpha = complete_space(n_orbitals, sector[1]);
    const auto beta = complete_space(n_orbitals, sector[2]);
    const int n_alpha = static_cast<int>(alpha.size());
    const int n_beta = static_cast<int>(beta.size());
    const int n_pairs = xmvb::vb::packed_active_pair_count(n_orbitals);
    xmvb::vb::OrthogonalActiveIntegrals integrals;
    integrals.one_electron = Eigen::MatrixXd::Identity(n_orbitals, n_orbitals);
    Eigen::MatrixXd factors(7, n_pairs);
    for (int j = 0; j < n_pairs; ++j)
      for (int i = 0; i < factors.rows(); ++i)
        factors(i, j) = 0.07 * std::sin(0.17 * (i + 1) * (j + 2));
    integrals.pair_kernel = factors.transpose() * factors;
    const xmvb::vb::DirectCiSigmaAction action(alpha, beta, integrals);
    for (int width : {1, 4}) {
      Eigen::MatrixXd coefficients(n_alpha, width * n_beta);
      for (int j = 0; j < coefficients.cols(); ++j)
        for (int i = 0; i < coefficients.rows(); ++i)
          coefficients(i, j) = std::sin(0.031 * (i + 1) * (j + 2));
      double checksum = action.apply(coefficients).squaredNorm(); // warmup
      constexpr int repeats = 5;
      auto start = Clock::now();
      for (int repeat = 0; repeat < repeats; ++repeat)
        checksum += action.apply(coefficients).squaredNorm();
      const double wall = std::chrono::duration<double>(Clock::now() - start).count();
      std::cout << "sigma," << n_orbitals << ',' << n_alpha << ',' << n_beta << ',' << width
                << ',' << repeats << ',' << wall << ',' << checksum << '\n';
      start = Clock::now();
      const auto adjoint = action.integral_adjoint(coefficients, coefficients);
      const double adjoint_wall = std::chrono::duration<double>(Clock::now() - start).count();
      std::cout << "adjoint," << n_orbitals << ',' << n_alpha << ',' << n_beta << ',' << width
                << ",1," << adjoint_wall << ','
                << adjoint.one_electron.squaredNorm() + adjoint.pair_kernel.squaredNorm()
                << '\n';
      Eigen::MatrixXd state = coefficients.array().cos().matrix();
      check_fused_response(action, coefficients, state);
      for (bool fused : {false, true}) {
        double total = 0.0;
        start = Clock::now();
        for (int repeat = 0; repeat < repeats; ++repeat) {
          if (fused) {
            const auto result = action.apply_with_adjoint(coefficients, state);
            total += result.sigma.squaredNorm() + result.integrals.one_electron.squaredNorm() +
                result.integrals.pair_kernel.squaredNorm();
          } else {
            const auto value = action.apply(coefficients);
            const auto gradient = action.integral_adjoint(coefficients, state);
            total += value.squaredNorm() + gradient.one_electron.squaredNorm() +
                gradient.pair_kernel.squaredNorm();
          }
        }
        const double elapsed = std::chrono::duration<double>(Clock::now() - start).count();
        std::cout << (fused ? "response_fused," : "response_separate,")
                  << n_orbitals << ',' << n_alpha << ',' << n_beta << ',' << width << ',' << repeats
                  << ',' << elapsed << ',' << total << '\n';
      }
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--benchmark-response") {
      benchmark_structure_response();
      return 0;
    }
    if (argc == 2 && std::string(argv[1]) == "--benchmark") {
      benchmark_sigma();
      return 0;
    }
    check_exterior_transform(7);
    check_exterior_transform(259);
    check_planner();
    check_sigma_action();
    check_sigma_tile_boundary();
    check_distinct_spin_sigma_action();
    check_fused_spin_sectors();
    check_demand_driven_structure_diagonal();
    std::cout << "orthogonal direct-CI planner and exterior transform: passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "orthogonal direct-CI test failed: " << error.what() << '\n';
    return 1;
  }
}
