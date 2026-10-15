#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>

#include "core/metric_solve.hpp"
#include "vbscf/determinants/pairs/same_spin_cache.hpp"
#include "vbscf/integrals/active/two_electron/construction/kernel.hpp"
#include "vbscf/structures/assembly/action.hpp"
#include "vbscf/structures/orthogonal_ci/restricted_metric.hpp"
#include "vbscf/structures/overlap/representation.hpp"

namespace {

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

void combinations(int n, int k, int first, std::vector<int>* current,
                  std::vector<std::vector<int>>* result) {
  if (static_cast<int>(current->size()) == k) {
    result->push_back(*current);
    return;
  }
  const int remaining = k - static_cast<int>(current->size());
  for (int orbital = first; orbital <= n - remaining; ++orbital) {
    current->push_back(orbital);
    combinations(n, k, orbital + 1, current, result);
    current->pop_back();
  }
}

std::vector<std::vector<int>> strings(int n, int k) {
  std::vector<std::vector<int>> result;
  std::vector<int> current;
  combinations(n, k, 0, &current, &result);
  return result;
}

struct Fixture {
  std::vector<std::vector<int>> alpha, beta;
  Eigen::MatrixXd upper;
  Eigen::MatrixXd selected;
  std::vector<std::vector<int>> product_alpha, product_beta;
  std::vector<std::vector<xmvb::vb::StructureExpansionTerm>> terms;
  xmvb::vb::SameSpinPairCacheContext topology;
  xmvb::vb::ActiveSpaceTwoElectronResult two_electron;
  Eigen::MatrixXd one_electron;
  std::vector<double> overlap;

  Fixture(int n, int n_alpha, int n_beta, Eigen::MatrixXd columns)
      : alpha(strings(n, n_alpha)), beta(strings(n, n_beta)),
        upper(Eigen::MatrixXd::Identity(n, n)), selected(std::move(columns)),
        one_electron(Eigen::MatrixXd::Zero(n, n)) {
    for (int row = 0; row < n; ++row) {
      upper(row, row) = 0.9 + 0.08 * row;
      for (int col = row + 1; col < n; ++col)
        upper(row, col) = 0.04 * (row + 1) * (col + 1);
    }
    const Eigen::MatrixXd s = upper.transpose() * upper;
    overlap.assign(s.data(), s.data() + s.size());
    terms.resize(alpha.size() * beta.size());
    for (const auto& a : alpha)
      for (const auto& b : beta) {
        const int product = static_cast<int>(product_alpha.size());
        product_alpha.push_back(a);
        product_beta.push_back(b);
        for (int col = 0; col < selected.cols(); ++col)
          if (selected(product, col) != 0.0)
            terms[product].push_back({col, selected(product, col)});
      }
    topology = xmvb::vb::build_same_spin_pair_topology(
        product_alpha, product_beta, n);
    two_electron.packed_active_two_electron_integrals.assign(
        xmvb::vb::packed_active_two_electron_integral_count(n), 0.0);
  }

  xmvb::vb::StructureAction action() const {
    return xmvb::vb::StructureAction(
        terms, selected.cols(), topology, overlap, one_electron, two_electron,
        upper.rows(), nullptr, xmvb::vb::StructureActionMode::ActionOnly);
  }
};

void compare_dense(Fixture& fixture, int expected_q) {
  auto action = fixture.action();
  require(action.storage().orthogonal_direct_ci,
          "test carrier did not select exact direct CI");
  auto inverse = xmvb::vb::RestrictedMetricInverse::make(
      action, fixture.alpha, fixture.beta, 64ULL * 1024ULL * 1024ULL);
  require(inverse != nullptr, "valid restricted spin sector was rejected");
  require(inverse->n_complement() == expected_q,
          "restricted metric chose the wrong complement dimension");
  auto prepared = inverse->prepare(fixture.upper);
  require(prepared != nullptr, "restricted metric preconditioner failed");
  xmvb::vb::RepresentationOverlap representation(
      action, fixture.alpha, fixture.beta, fixture.upper.rows());
  const auto prepared_metric = representation.prepare_metric(fixture.upper);
  const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(
      fixture.selected.cols(), fixture.selected.cols());
  const Eigen::MatrixXd metric = representation.apply(
      fixture.upper, fixture.upper, identity);
  Eigen::MatrixXd rhs(fixture.selected.cols(), 11);
  for (int col = 0; col < rhs.cols(); ++col)
    for (int row = 0; row < rhs.rows(); ++row)
      rhs(row, col) = std::sin(0.21 * (row + 1) * (col + 2));
  const Eigen::MatrixXd preconditioned = prepared->precondition(rhs);
  require((preconditioned.leftCols(1) -
           prepared->precondition(rhs.leftCols(1))).norm() < 1.0e-12,
          "cached principal inverse changed across repeated RHS");
  require(rhs.col(0).dot(preconditioned.col(0)) > 0.0,
          "principal inverse is not a positive preconditioner");
  Eigen::MatrixXd solution(rhs.rows(), rhs.cols());
  for (int col = 0; col < rhs.cols(); ++col) {
    xmvb::core::MetricSolveStats stats;
    solution.col(col) = xmvb::core::solve_metric(
        Eigen::VectorXd(rhs.col(col)),
        [&](const Eigen::VectorXd& x) {
          return Eigen::VectorXd(prepared_metric(x).col(0));
        },
        [&](const Eigen::VectorXd& x) {
          return Eigen::VectorXd(prepared->precondition(x).col(0));
        }, &stats);
    require(stats.certifications > 0 &&
                stats.relative_residual <= std::sqrt(std::numeric_limits<double>::epsilon()),
            "restricted PCG did not certify the selected metric residual");
  }
  const Eigen::MatrixXd reference = metric.fullPivLu().solve(rhs);
  require((solution - reference).norm() <
              2.0e-10 * std::max(1.0, reference.norm()),
          "restricted metric differs from independently solved selected Gram");
  require((metric * solution - rhs).norm() <
              2.0e-10 * std::max(1.0, rhs.norm()),
          "restricted metric failed selected-space residual");
}

void singlet_case() {
  Eigen::MatrixXd selected = Eigen::MatrixXd::Zero(9, 4);
  selected(0, 0) = 1.0;
  selected(4, 1) = 1.0;
  selected(8, 2) = 1.0;
  selected(1, 3) = selected(3, 3) = 1.0;
  Fixture fixture(3, 1, 1, selected);
  compare_dense(fixture, 2);

  Eigen::MatrixXd full = Eigen::MatrixXd::Zero(9, 6);
  full.leftCols(4) = selected;
  full(2, 4) = full(6, 4) = 1.0;
  full(5, 5) = full(7, 5) = 1.0;
  Eigen::MatrixXd compound(9, 9);
  for (int target_a = 0; target_a < 3; ++target_a)
    for (int target_b = 0; target_b < 3; ++target_b)
      for (int source_a = 0; source_a < 3; ++source_a)
        for (int source_b = 0; source_b < 3; ++source_b)
          compound(target_a * 3 + target_b, source_a * 3 + source_b) =
              fixture.upper(target_a, source_a) *
              fixture.upper(target_b, source_b);
  const Eigen::MatrixXd carrier_metric = compound.transpose() * compound;
  const Eigen::MatrixXd selected_metric =
      selected.transpose() * carrier_metric * selected;
  const Eigen::MatrixXd full_metric = full.transpose() * carrier_metric * full;
  const Eigen::MatrixXd naive_pp =
      full_metric.inverse().topLeftCorner(4, 4);
  auto action = fixture.action();
  auto inverse = xmvb::vb::RestrictedMetricInverse::make(
      action, fixture.alpha, fixture.beta, 64ULL * 1024ULL * 1024ULL);
  require(inverse != nullptr, "singlet preconditioner was not built");
  auto prepared = inverse->prepare(fixture.upper);
  require((prepared->precondition(Eigen::Matrix4d::Identity()) - naive_pp).norm() < 1.0e-10,
          "principal full inverse differs from explicit N_PP");
  require((selected_metric.inverse() - naive_pp).norm() > 1.0e-9,
          "test failed to distinguish the restricted inverse from N_PP");

  selected.col(3) *= -2.0;
  Fixture scaled(3, 1, 1, selected);
  scaled.upper(0, 0) = 0.35;
  scaled.upper(2, 2) = 2.8;
  compare_dense(scaled, 2);

  selected(3, 3) = 0.0;  // Not annihilated by S_+.
  Fixture impure(3, 1, 1, selected);
  auto impure_action = impure.action();
  require(xmvb::vb::RestrictedMetricInverse::make(
              impure_action, impure.alpha, impure.beta, 64ULL * 1024ULL * 1024ULL) ==
              nullptr,
          "spin-impure selected column was admitted");

  selected(3, 3) = -2.0;
  selected(0, 3) = 0.2;  // One column now crosses occupation blocks.
  Fixture nonlocal(3, 1, 1, selected);
  auto nonlocal_action = nonlocal.action();
  require(xmvb::vb::RestrictedMetricInverse::make(
              nonlocal_action, nonlocal.alpha, nonlocal.beta,
              64ULL * 1024ULL * 1024ULL) == nullptr,
          "occupation-mixed selected column was admitted");

  auto valid_action = fixture.action();
  require(xmvb::vb::RestrictedMetricInverse::make(
              valid_action, fixture.alpha, fixture.beta, 1024) == nullptr,
          "restricted metric ignored its fixed byte budget");
}

void doublet_case() {
  // Three orbitals, two alpha and one beta: six one-doubly-occupied
  // highest-weight products and two three-open-shell doublets.
  Eigen::MatrixXd selected = Eigen::MatrixXd::Zero(9, 7);
  int col = 0;
  const auto alpha = strings(3, 2);
  const auto beta = strings(3, 1);
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      if (std::find(alpha[a].begin(), alpha[a].end(), beta[b][0]) != alpha[a].end())
        selected(a * 3 + b, col++) = 1.0;
  require(col == 6, "doublet fixture has wrong paired-product count");
  // S_+ signs for (01;2) and (02;1) are + and -, respectively.
  selected(0 * 3 + 2, 6) = 1.0;
  selected(1 * 3 + 1, 6) = 1.0;
  Fixture fixture(3, 2, 1, selected);
  compare_dense(fixture, 1);
}

}  // namespace

int main() {
  singlet_case();
  doublet_case();
  std::cout << "restricted metric checks passed\n";
}
