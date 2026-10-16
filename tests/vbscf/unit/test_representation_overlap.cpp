#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/QR>

#include "vbscf/determinants/pairs/same_spin_cache.hpp"
#include "vbscf/integrals/active/two_electron/construction/kernel.hpp"
#include "vbscf/structures/assembly/action.hpp"
#include "vbscf/structures/overlap/representation.hpp"
#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"

namespace {

using namespace xmvb::vb;

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

double relative_error(const Eigen::MatrixXd& actual,
                      const Eigen::MatrixXd& reference) {
  return (actual - reference).norm() /
         std::max(1.0, reference.norm());
}

double dot(const Eigen::MatrixXd& left, const Eigen::MatrixXd& right) {
  return (left.array() * right.array()).sum();
}

Eigen::MatrixXd upper(double diagonal_shift, double off_diagonal_scale) {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(4, 4);
  for (int i = 0; i < 4; ++i) {
    result(i, i) = 1.0 + diagonal_shift + 0.09 * i;
    for (int j = i + 1; j < 4; ++j) {
      result(i, j) = off_diagonal_scale *
                     std::sin(0.61 * (i + 1) * (j + 2));
    }
  }
  return result;
}

Eigen::MatrixXd upper_direction(double scale) {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(4, 4);
  for (int i = 0; i < 4; ++i) {
    for (int j = i; j < 4; ++j) {
      result(i, j) = scale * std::cos(0.37 * (i + 2) * (j + 1));
    }
  }
  return result;
}

Eigen::MatrixXd selected_expansion(int n_alpha, int n_beta) {
  Eigen::MatrixXd w(n_alpha * n_beta, 3);
  for (int product = 0; product < w.rows(); ++product) {
    w(product, 0) = std::sin(0.41 * (product + 1));
    w(product, 1) = std::cos(0.29 * (product + 2));
    w(product, 2) = std::sin(0.17 * (product + 3)) +
                    0.37 * std::cos(0.53 * (product + 1));
  }
  require(Eigen::FullPivLU<Eigen::MatrixXd>(w).rank() == 3,
          "selected expansion lost column rank");
  return w;
}

StructureAction make_action(
    const std::vector<std::vector<int>>& alpha,
    const std::vector<std::vector<int>>& beta,
    const Eigen::MatrixXd& w) {
  const int n_beta = static_cast<int>(beta.size());
  int n_orbitals=0;
  for (const auto& strings : {alpha,beta})
    for (const auto& occupied : strings)
      for (int p:occupied) n_orbitals=std::max(n_orbitals,p+1);
  std::vector<std::vector<int>> determinant_alpha, determinant_beta;
  std::vector<std::vector<StructureExpansionTerm>> terms(w.rows());
  for (int a = 0; a < static_cast<int>(alpha.size()); ++a) {
    for (int b = 0; b < n_beta; ++b) {
      const int product = a * n_beta + b;
      determinant_alpha.push_back(alpha[a]);
      determinant_beta.push_back(beta[b]);
      for (int structure = 0; structure < 3; ++structure) {
        if (w(product, structure) != 0.0)
          terms[product].push_back({structure, w(product, structure)});
      }
    }
  }
  const auto topology = build_same_spin_pair_topology(
      determinant_alpha, determinant_beta, n_orbitals);
  const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(n_orbitals,n_orbitals);
  const std::vector<double> packed_overlap(identity.data(),
                                            identity.data() + identity.size());
  ActiveSpaceTwoElectronResult two_electron;
  two_electron.packed_active_two_electron_integrals.assign(
      packed_active_two_electron_integral_count(n_orbitals), 0.0);
  return StructureAction(
      terms, 3, topology, packed_overlap, Eigen::MatrixXd::Zero(n_orbitals,n_orbitals),
      two_electron, n_orbitals, nullptr, StructureActionMode::ActionOnly);
}

// Independent small-system oracle: enumerate each occupied minor directly.
Eigen::MatrixXd compound(
    const Eigen::MatrixXd& kernel,
    const std::vector<std::vector<int>>& occupied) {
  const int n = static_cast<int>(occupied.size());
  const int electrons = static_cast<int>(occupied.front().size());
  Eigen::MatrixXd result(n, n);
  for (int left = 0; left < n; ++left) {
    for (int right = 0; right < n; ++right) {
      Eigen::MatrixXd minor(electrons, electrons);
      for (int row = 0; row < electrons; ++row) {
        for (int col = 0; col < electrons; ++col) {
          minor(row, col) =
              kernel(occupied[left][row], occupied[right][col]);
        }
      }
      result(left, right) = minor.determinant();
    }
  }
  return result;
}

Eigen::MatrixXd product_compound(
    const Eigen::MatrixXd& kernel,
    const std::vector<std::vector<int>>& alpha,
    const std::vector<std::vector<int>>& beta) {
  const Eigen::MatrixXd a = compound(kernel, alpha);
  const Eigen::MatrixXd b = compound(kernel, beta);
  Eigen::MatrixXd product(a.rows() * b.rows(), a.cols() * b.cols());
  for (int al = 0; al < a.rows(); ++al) {
    for (int bl = 0; bl < b.rows(); ++bl) {
      for (int ar = 0; ar < a.cols(); ++ar) {
        for (int br = 0; br < b.cols(); ++br) {
          product(al * b.rows() + bl, ar * b.cols() + br) =
              a(al, ar) * b(bl, br);
        }
      }
    }
  }
  return product;
}

Eigen::MatrixXd reference_cross(
    const Eigen::MatrixXd& left,
    const Eigen::MatrixXd& right,
    const Eigen::MatrixXd& w,
    const std::vector<std::vector<int>>& alpha,
    const std::vector<std::vector<int>>& beta) {
  return w.transpose() * product_compound(left.transpose() * right,
                                          alpha, beta) * w;
}

void check_case(const std::vector<std::vector<int>>& alpha,
                const std::vector<std::vector<int>>& beta,
                bool full_ambient, bool sparse = false, bool matching = false) {
  Eigen::MatrixXd w = selected_expansion(alpha.size(), beta.size());
  if (sparse) {
    w.setZero();
    // Arbitrary three structures with strictly restricted carrier support.
    // Never seed the highest product, which would make the support full.
    w(0, 0) = 1.0;
    w(1, 1) = -0.7;
    w(beta.size(), 1) = 0.2;
    w(beta.size() + 1, 2) = 0.9;
  }
  if (matching) {
    w.setZero();
    for (int b=0;b<static_cast<int>(beta.size());++b)
      for (int s=0;s<3;++s)
        w((b%alpha.size())*beta.size()+b,s)=std::sin(0.37*(b+1)*(s+2));
  }
  const auto action = make_action(alpha, beta, w);
  const RepresentationOverlap overlap(action, alpha, beta, 4);
  const Eigen::MatrixXd left = upper(0.04, 0.29);
  const Eigen::MatrixXd right = upper(-0.08, -0.21);
  const Eigen::MatrixXd dl = upper_direction(0.46);
  const Eigen::MatrixXd dr = upper_direction(-0.33);
  const Eigen::MatrixXd kernel = left.transpose() * right;
  require((kernel - kernel.transpose()).norm() > 0.1,
          "cross-kernel fixture is accidentally symmetric");
  Eigen::Vector3d x(0.6, -0.4, 0.8), y(-0.3, 0.7, 0.5);
  Eigen::MatrixXd block(3, 2);
  block.col(0) = y;
  block.col(1) = x;

  const Eigen::MatrixXd reference =
      reference_cross(left, right, w, alpha, beta);
  require(relative_error(overlap.apply(left, right, block),
                         reference * block) < 2e-12,
          "cross-Gram block action differs from explicit occupied minors");
  const double scalar = x.dot(overlap.apply(left, right, y).col(0));
  require(std::abs(scalar - x.dot(reference * y)) < 2e-12,
          "cross-Gram scalar differs from explicit occupied minors");

  const Eigen::MatrixXd metric_reference =
      reference_cross(right, right, w, alpha, beta);
  const auto prepared = overlap.prepare_metric(right);
  require(relative_error(prepared(block), metric_reference * block) < 2e-12 &&
              relative_error(prepared(block),
                             overlap.apply(right, right, block)) < 2e-12,
          "prepared metric multi-RHS action differs from cross-Gram action");
  require(relative_error(prepared(Eigen::Matrix3d::Identity()),
                         metric_reference) < 2e-12,
          "prepared metric matrix differs from explicit occupied minors");
  require(relative_error(prepared(x), metric_reference * x) < 2e-12,
          "prepared metric vector action differs from explicit occupied minors");
  require(Eigen::LLT<Eigen::MatrixXd>(metric_reference).info() == Eigen::Success &&
              x.dot(prepared(x).col(0)) > 0.0,
          "prepared representation metric is not positive definite");
  Eigen::MatrixXd mutable_chart = right;
  const auto frozen = overlap.prepare_metric(mutable_chart);
  mutable_chart(0, 0) += 0.2;
  require(relative_error(frozen(block), metric_reference * block) < 2e-12,
          "prepared metric did not retain its fixed chart");

  const Eigen::MatrixXd reverse =
      overlap.apply(right, left, Eigen::Matrix3d::Identity());
  require(relative_error(reverse, reference.transpose()) < 2e-12,
          "left/right cross-Gram transpose symmetry failed");
  require(std::abs(scalar - y.dot(overlap.apply(right, left, x).col(0))) < 2e-12,
          "left/right scalar cross symmetry failed");

  constexpr double h = 1e-6;
  const Eigen::MatrixXd finite_action_direction =
      (overlap.apply(left + h * dl, right + h * dr, block) -
       overlap.apply(left - h * dl, right - h * dr, block)) / (2.0 * h);
  const Eigen::MatrixXd analytic_action_direction =
      overlap.apply_direction(left, right, dl, dr, block);
  require(relative_error(analytic_action_direction,
                         finite_action_direction) < 2e-8,
          "cross-Gram directional block action fails finite differences");

  if (full_ambient) {
    ExteriorOrbitalTransform la(alpha, left), lb(beta, left);
    ExteriorOrbitalTransform ra(alpha, right), rb(beta, right);
    const auto dla = la.direction(dl), dlb = lb.direction(dl);
    const auto dra = ra.direction(dr), drb = rb.direction(dr);
    const std::size_t one_column = 8 * sizeof(double) * (alpha.size() + beta.size());
    require(relative_error(action.cross_overlap(block, la, lb, ra, rb,
        nullptr, nullptr, nullptr, nullptr, one_column), reference * block) < 2e-12,
        "bounded cross-Gram panels differ from independent minors");
    require(relative_error(action.cross_overlap(block, la, lb, ra, rb,
        &dla, &dlb, &dra, &drb, one_column), analytic_action_direction) < 2e-12,
        "one-column derivative panels differ from wide panels");
    const auto tiled = action.cross_overlap_adjoint(x, y, la, lb, ra, rb,
        nullptr, nullptr, 2 * one_column);
    const auto wide = overlap.adjoint(left, right, x, y);
    require(relative_error(tiled.left, wide.left) < 2e-12 &&
        relative_error(tiled.right, wide.right) < 2e-12,
        "one-column adjoint differs from wide panels");
    const auto dtiled = action.cross_overlap_adjoint(x, y, la, lb, ra, rb,
        &dl, &dr, 2 * one_column);
    const auto dwide = overlap.adjoint_direction(left, right, dl, dr, x, y);
    require(relative_error(dtiled.left, dwide.left) < 2e-12 &&
        relative_error(dtiled.right, dwide.right) < 2e-12,
        "one-column adjoint derivative differs from wide panels");
  }

  const UpperPair adjoint = overlap.adjoint(left, right, x, y);
  const double adjoint_dot = dot(adjoint.left, dl) + dot(adjoint.right, dr);
  const double scalar_difference =
      (x.dot(overlap.apply(left + h * dl, right + h * dr, y).col(0)) -
       x.dot(overlap.apply(left - h * dl, right - h * dr, y).col(0))) /
      (2.0 * h);
  require(std::abs(adjoint_dot - scalar_difference) <
              2e-8 * std::max(1.0, std::abs(scalar_difference)),
          "cross-Gram upper adjoint fails scalar directional difference");
  require(std::abs(adjoint_dot -
                   x.dot(overlap.apply_direction(left, right, dl, dr, y).col(0))) <
              2e-10,
          "cross-Gram action and adjoint violate dot-product duality");
  const UpperPair swapped = overlap.adjoint(right, left, y, x);
  require(relative_error(adjoint.left, swapped.right) < 2e-12 &&
              relative_error(adjoint.right, swapped.left) < 2e-12,
          "left/right adjoints violate cross symmetry");

  const UpperPair analytic_adjoint_direction =
      overlap.adjoint_direction(left, right, dl, dr, x, y);
  const UpperPair plus = overlap.adjoint(left + h * dl, right + h * dr, x, y);
  const UpperPair minus = overlap.adjoint(left - h * dl, right - h * dr, x, y);
  require(relative_error(analytic_adjoint_direction.left,
                         (plus.left - minus.left) / (2.0 * h)) < 2e-7 &&
              relative_error(analytic_adjoint_direction.right,
                             (plus.right - minus.right) / (2.0 * h)) < 2e-7,
          "cross-Gram adjoint direction fails finite differences");

  // At identity, many pair minors vanish exactly. Cofactor derivatives of
  // those singular minors must still contribute to the chart adjoint.
  const Eigen::Matrix4d identity = Eigen::Matrix4d::Identity();
  require(std::abs(compound(identity, alpha)(0, 1)) < 1e-15,
          "identity fixture does not contain a singular pair minor");
  const UpperPair identity_adjoint = overlap.adjoint(identity, identity, x, y);
  const double identity_difference =
      (x.dot(overlap.apply(identity + h * dl, identity, y).col(0)) -
       x.dot(overlap.apply(identity - h * dl, identity, y).col(0))) /
      (2.0 * h);
  require(std::abs(dot(identity_adjoint.left, dl) - identity_difference) < 2e-8,
          "singular occupied minors lost their derivative contribution");

  if (full_ambient) {
    // The carrier is complete, but the three selected columns do not close
    // under a finite orbital change.
    const Eigen::MatrixXd f_left =
        product_compound(left, alpha, beta) * w;
    const Eigen::MatrixXd f_right =
        product_compound(right, alpha, beta) * w;
    const Eigen::VectorXd target = f_right * y;
    const Eigen::VectorXd residual =
        target - f_left * f_left.colPivHouseholderQr().solve(target);
    require(residual.norm() > 1e-3,
            "full-ambient selected-space fixture accidentally closes");
  }
}

/** @brief Both exact contraction orders, including one-sided chart derivatives. */
void check_broad_moment_budget() {
  std::vector<std::vector<int>> spin;
  for (int p=0;p<5;++p) for (int q=p+1;q<5;++q) spin.push_back({p,q});
  const auto w=selected_expansion(spin.size(),spin.size());
  const auto action=make_action(spin,spin,w);
  Eigen::MatrixXd l=Eigen::MatrixXd::Identity(5,5),r=1.2*l;
  Eigen::MatrixXd dl=Eigen::MatrixXd::Zero(5,5),dr=dl;
  for (int q=0;q<5;++q) for (int p=0;p<=q;++p) {
    if (p<q) { l(p,q)=0.07*std::sin(p+2*q+1); r(p,q)=0.08*std::cos(2*p+q+1); }
    dl(p,q)=0.03*std::cos(p+q+2); dr(p,q)=0.04*std::sin(p+2*q+2);
  }
  const ExteriorOrbitalTransform la(spin,l),ra(spin,r);
  const Eigen::Vector3d x(0.4,-0.7,0.8),y(-0.2,0.6,0.3);
  const auto primal=action.cross_overlap_adjoint(x,y,la,la,ra,ra);
  const auto factor=action.cross_overlap_adjoint(x,y,la,la,ra,ra,nullptr,nullptr,1280);
  require(relative_error(primal.left,factor.left)<2e-12 &&
              relative_error(primal.right,factor.right)<2e-12,
          "broad one-body moment differs from the bounded source adjoint");
  const Eigen::MatrixXd zero=Eigen::MatrixXd::Zero(5,5);
  for (bool left_direction : {false,true}) for (bool right_direction : {false,true}) {
    const auto& a=left_direction ? dl:zero;
    const auto& b=right_direction ? dr:zero;
    const auto dense=action.cross_overlap_adjoint(x,y,la,la,ra,ra,&a,&b);
    const auto narrow=action.cross_overlap_adjoint(x,y,la,la,ra,ra,&a,&b,2560);
    require(relative_error(dense.left,narrow.left)<3e-12 &&
                relative_error(dense.right,narrow.right)<3e-12,
            "broad moment derivative differs from the bounded source derivative");
    const auto da=la.direction(a),db=ra.direction(b);
    const auto image=action.cross_overlap(y,la,la,ra,ra,&da,&da,&db,&db);
    const auto tiled=action.cross_overlap(y,la,la,ra,ra,&da,&da,&db,&db,1280);
    require(relative_error(image,tiled)<3e-12,
            "broad cross-overlap jet differs from bounded source panels");
  }
}

}  // namespace

int main() {
  try {
    check_broad_moment_budget();
    const std::vector<std::vector<int>> full_alpha = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    const std::vector<std::vector<int>> full_beta = {
        {0}, {1}, {2}, {3}};
    check_case(full_alpha, full_beta, true);
    check_case(full_beta, full_beta, true);
    check_case(full_alpha, full_beta, true, true);
    check_case(full_beta, full_beta, true, true);
    check_case(full_alpha, full_beta, true, false, true);
    check_case(full_beta, full_beta, true, false, true);
    check_case({{0, 1}, {0, 2}, {1, 3}}, {{0}, {2}, {3}}, false);
    std::cout << "Selected-space representation cross-Gram and derivatives: passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
