#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "vbscf/optimization/neo/ciah.hpp"

namespace {

using xmvb::vb::CiahOptions;
using xmvb::vb::CiahResult;
using xmvb::vb::CiahWorkspace;
using xmvb::vb::NeoProblem;

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

void require_close(
    const Eigen::VectorXd& actual, const Eigen::VectorXd& expected,
    double tolerance, const std::string& message) {
  require((actual - expected).stableNorm() <=
              tolerance * std::max(1.0, expected.stableNorm()),
          message);
}

double dual_norm(
    const Eigen::MatrixXd& metric, const Eigen::VectorXd& covector) {
  return std::sqrt(covector.dot(metric.ldlt().solve(covector)));
}

NeoProblem dense_problem(
    const Eigen::MatrixXd& hessian, const Eigen::MatrixXd& metric,
    const Eigen::VectorXd& gradient, int* actions = nullptr) {
  return NeoProblem(
      gradient,
      [hessian, actions](const Eigen::VectorXd& vector) {
        if (actions) ++*actions;
        return Eigen::VectorXd(hessian * vector);
      },
      [metric](const Eigen::VectorXd& vector) {
        return Eigen::VectorXd(metric * vector);
      },
      [metric](const Eigen::VectorXd& covector) {
        return dual_norm(metric, covector);
      });
}

void check_images_and_residual(
    const Eigen::MatrixXd& hessian, const Eigen::MatrixXd& metric,
    const Eigen::VectorXd& gradient, const CiahResult& result) {
  require_close(result.hessian_step, hessian * result.step, 2.0e-12,
                "CIAH returned an inconsistent Hessian image");
  require_close(result.metric_step, metric * result.step, 2.0e-12,
                "CIAH returned an inconsistent metric image");
  const Eigen::VectorXd expected = result.reference_weight *
      (gradient + hessian * result.step -
       result.eigenvalue * metric * result.step);
  // The eigenvector phase is free, so the stored residual may have either sign.
  const double error = std::min(
      (result.residual - expected).stableNorm(),
      (result.residual + expected).stableNorm());
  require(error <= 2.0e-12 * std::max(1.0, expected.stableNorm()),
          "CIAH returned an inconsistent full-space AH residual");
  require(std::abs(result.residual_norm -
                   dual_norm(metric, result.residual)) <= 2.0e-12,
          "CIAH residual norm is not the physical dual norm");
}

void check_full_space_root(
    const Eigen::MatrixXd& hessian, const Eigen::MatrixXd& metric,
    const Eigen::VectorXd& gradient, const CiahResult& result) {
  const Eigen::Index n = gradient.size();
  Eigen::MatrixXd augmented_hessian = Eigen::MatrixXd::Zero(n + 1, n + 1);
  augmented_hessian.block(0, 1, 1, n) = gradient.transpose();
  augmented_hessian.block(1, 0, n, 1) = gradient;
  augmented_hessian.bottomRightCorner(n, n) = hessian;
  Eigen::MatrixXd augmented_metric =
      Eigen::MatrixXd::Identity(n + 1, n + 1);
  augmented_metric.bottomRightCorner(n, n) = metric;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(
      augmented_hessian, augmented_metric);
  require(spectrum.info() == Eigen::Success,
          "dense generalized AH reference failed");
  Eigen::Index selected = -1;
  for (Eigen::Index i = 0; i < n + 1; ++i) {
    if (std::abs(spectrum.eigenvectors()(0, i)) > 0.1) {
      selected = i;
      break;
    }
  }
  require(selected >= 0, "dense generalized AH has no admissible root");
  const double alpha = spectrum.eigenvectors()(0, selected);
  const Eigen::VectorXd expected_step =
      spectrum.eigenvectors().col(selected).tail(n) / alpha;
  require(std::abs(result.eigenvalue - spectrum.eigenvalues()[selected]) <
              2.0e-11 * std::max(1.0, std::abs(result.eigenvalue)),
          "CIAH selected the wrong generalized AH root");
  require_close(result.step, expected_step, 2.0e-10,
                "CIAH step differs from dense generalized AH reference");
}

void check_anisotropic_spd() {
  Eigen::Matrix3d hessian;
  hessian << 4.0, 0.3, -0.2,
             0.3, 2.5, 0.4,
            -0.2, 0.4, 5.0;
  Eigen::Matrix3d metric;
  metric << 1.0e-3, 0.0, 0.0,
            0.0, 7.0, 0.2,
            0.0, 0.2, 2.0e2;
  const Eigen::Vector3d gradient(0.3, -1.2, 0.7);
  int actions = 0;
  const NeoProblem problem = dense_problem(hessian, metric, gradient, &actions);
  CiahWorkspace workspace(problem, Eigen::VectorXd::Ones(3));
  CiahResult result;
  do {
    result = workspace.next(gradient);
    require(result.new_hvp_actions <= 1,
            "CIAH performed more than one Hessian action in a call");
    require(actions == result.basis_size,
            "CIAH repeated a cached Hessian action");
    check_images_and_residual(hessian, metric, gradient, result);
  } while (!result.converged &&
           (!result.exhaustion || result.new_hvp_actions != 0));
  if (!result.converged) result = workspace.next(gradient);
  require(result.converged,
          "anisotropic generalized CIAH failed to converge in full space");
  check_full_space_root(hessian, metric, gradient, result);
}

void check_indefinite_root_selection() {
  Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();
  hessian.diagonal() << -100.0, 2.0, 4.0;
  const Eigen::Matrix3d metric = Eigen::Matrix3d::Identity();
  const Eigen::Vector3d gradient(0.0, 1.0, 0.0);
  const NeoProblem problem = dense_problem(hessian, metric, gradient);
  CiahWorkspace workspace(problem, Eigen::VectorXd::Unit(3, 0));
  CiahResult result;
  do {
    result = workspace.next(gradient);
  } while (!result.converged &&
           (!result.exhaustion || result.new_hvp_actions != 0));
  if (!result.converged) result = workspace.next(gradient);
  require(result.converged, "indefinite CIAH root did not converge");
  check_full_space_root(hessian, metric, gradient, result);
  require(result.eigenvalue > -1.0,
          "CIAH selected a Hessian-only root with zero reference weight");
  check_images_and_residual(hessian, metric, gradient, result);
}

void check_gradient_refresh_and_preconditioner() {
  Eigen::Matrix2d hessian;
  hessian << 3.0, 0.0,
             0.0, 6.0;
  const Eigen::Matrix2d metric = Eigen::Matrix2d::Identity();
  const Eigen::Vector2d old_gradient(1.0, 0.0);
  const Eigen::Vector2d new_gradient(0.2, 1.0);
  int actions = 0;
  int preconditioner_calls = 0;
  double last_eigenvalue = 0.0;
  CiahOptions options;
  options.residual_preconditioner =
      [&](const Eigen::VectorXd& residual, double eigenvalue) {
        ++preconditioner_calls;
        last_eigenvalue = eigenvalue;
        return Eigen::VectorXd(residual);
      };
  const NeoProblem problem = dense_problem(
      hessian, metric, old_gradient, &actions);
  CiahWorkspace workspace(problem, Eigen::VectorXd::Unit(2, 0), options);
  const CiahResult first = workspace.next(old_gradient);
  require(first.residual_norm < 1.0e-12 && actions == 1,
          "one-dimensional AH residual should vanish with one action");
  const CiahResult settled = workspace.next(old_gradient);
  require(settled.converged && settled.new_hvp_actions == 0,
          "unchanged AH root did not satisfy eigenvalue-change convergence");
  const CiahResult refreshed = workspace.next(new_gradient);
  require(refreshed.new_hvp_actions == 0 && actions == 1 &&
              refreshed.basis_size == first.basis_size,
          "gradient refresh recomputed a cached Hessian image");
  require(!refreshed.converged,
          "gradient refresh unexpectedly remained converged");
  check_images_and_residual(hessian, metric, new_gradient, refreshed);
  const CiahResult expanded = workspace.next(new_gradient);
  require(expanded.new_hvp_actions == 1 && actions == 2 &&
              preconditioner_calls == 1,
          "residual expansion did not use exactly one shifted action");
  require(std::abs(last_eigenvalue - refreshed.eigenvalue) < 1.0e-12,
          "residual preconditioner received the wrong eigenvalue");
  const CiahResult converged = workspace.next(new_gradient);
  require(converged.converged,
          "refreshed full-space AH did not converge");
  check_full_space_root(hessian, metric, new_gradient, expanded);
  check_images_and_residual(hessian, metric, new_gradient, expanded);
  const CiahResult repeated = workspace.next(new_gradient);
  require(repeated.new_hvp_actions == 0 && actions == 2,
          "repeated CIAH query recomputed a Hessian image");
}

void check_metric_validation() {
  const Eigen::Matrix2d hessian = Eigen::Matrix2d::Identity();
  Eigen::Matrix2d metric = Eigen::Matrix2d::Zero();
  metric(0, 0) = 1.0;
  const Eigen::Vector2d gradient(1.0, 0.0);
  const NeoProblem problem = dense_problem(hessian, metric, gradient);
  CiahWorkspace workspace(problem, Eigen::VectorXd::Unit(2, 1));
  bool rejected = false;
  try {
    workspace.next(gradient);
  } catch (const std::runtime_error&) {
    rejected = true;
  }
  require(rejected, "CIAH accepted a rank-deficient projected metric");
}

}  // namespace

int main() {
  try {
    check_anisotropic_spd();
    check_indefinite_root_selection();
    check_gradient_refresh_and_preconditioner();
    check_metric_validation();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
