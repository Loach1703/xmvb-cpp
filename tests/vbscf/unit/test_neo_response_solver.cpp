#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "vbscf/optimization/neo/response_solver.hpp"
#include "vbscf/optimization/trust_region/spectral.hpp"

namespace {

using xmvb::vb::NeoOptions;
using xmvb::vb::NeoStopReason;
using xmvb::vb::ResponseNeoDirection;
using xmvb::vb::ResponseNeoProblem;
using xmvb::vb::ResponseNeoResult;
using xmvb::vb::ResponseNeoWorkspace;

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

void require_close(
    const Eigen::VectorXd& actual,
    const Eigen::VectorXd& reference,
    double tolerance,
    const std::string& message) {
  require((actual - reference).stableNorm() <=
              tolerance * std::max(1.0, reference.stableNorm()),
          message);
}

double physical_dual_norm(
    const Eigen::MatrixXd& metric,
    const Eigen::VectorXd& covector) {
  return std::sqrt(covector.dot(metric.ldlt().solve(covector)));
}

ResponseNeoProblem dense_problem(
    const Eigen::MatrixXd& a,
    const Eigen::MatrixXd& b,
    const Eigen::MatrixXd& c,
    const Eigen::MatrixXd& metric,
    const Eigen::VectorXd& gradient) {
  const Eigen::MatrixXd metric_inverse = metric.ldlt().solve(
      Eigen::MatrixXd::Identity(metric.rows(), metric.cols()));
  return ResponseNeoProblem(
      gradient,
      c.rows(),
      [a, b](const Eigen::VectorXd& p) {
        return ResponseNeoDirection{a * p, b * p};
      },
      [b, c](const Eigen::VectorXd& q) {
        return ResponseNeoDirection{b.transpose() * q, c * q};
      },
      [metric](const Eigen::VectorXd& p) { return metric * p; },
      [metric_inverse](const Eigen::VectorXd& covector) {
        return std::sqrt(covector.dot(metric_inverse * covector));
      });
}

Eigen::VectorXd explicit_orbital_step(
    const Eigen::MatrixXd& relaxed,
    const Eigen::MatrixXd& metric,
    const Eigen::VectorXd& gradient,
    double radius) {
  Eigen::LLT<Eigen::MatrixXd> factor(metric);
  require(factor.info() == Eigen::Success, "explicit metric factor failed");
  const Eigen::MatrixXd whitening = factor.matrixU().solve(
      Eigen::MatrixXd::Identity(metric.rows(), metric.cols()));
  const Eigen::MatrixXd white_hessian =
      whitening.transpose() * relaxed * whitening;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(white_hessian);
  require(spectrum.info() == Eigen::Success, "explicit spectrum failed");
  const auto solution = xmvb::vb::solve_spectral_trust_region(
      spectrum.eigenvalues(),
      spectrum.eigenvectors().transpose() *
          whitening.transpose() * gradient,
      radius);
  return whitening * spectrum.eigenvectors() * solution.step;
}

void verify_coupled_residual(
    const Eigen::MatrixXd& a,
    const Eigen::MatrixXd& b,
    const Eigen::MatrixXd& c,
    const Eigen::MatrixXd& metric,
    const Eigen::VectorXd& gradient,
    double radius,
    const ResponseNeoResult& result) {
  require(result.converged(), "response NEO did not converge");
  const Eigen::VectorXd orbital = gradient +
      a * result.step.orbital + b.transpose() * result.step.structure +
      result.shift * metric * result.step.orbital;
  const Eigen::VectorXd structure =
      b * result.step.orbital + c * result.step.structure;
  const double residual = std::hypot(
      physical_dual_norm(metric, orbital), structure.stableNorm());
  const double scale = std::max({
      1.0, physical_dual_norm(metric, gradient),
      physical_dual_norm(metric, a * result.step.orbital),
      physical_dual_norm(metric, b.transpose() * result.step.structure)});
  require(residual <= 2.0e-9 * scale,
          "full coupled KKT residual is not converged");
  require(std::abs(result.residual_norm - residual) <=
              2.0e-9 * std::max(1.0, residual),
          "returned coupled residual has the wrong physical norm");
  require(structure.stableNorm() <= 2.0e-9 * scale,
          "structure stationarity is not converged");
  require(std::sqrt(result.step.orbital.dot(
              metric * result.step.orbital)) <=
              radius * (1.0 + 2.0e-10),
          "orbital step violates its trust region");
}

void check_spd_structure_and_explicit_schur() {
  Eigen::Matrix3d a;
  a << 3.0, 0.4, -0.2,
       0.4, 2.2, 0.3,
      -0.2, 0.3, 1.7;
  Eigen::Matrix<double, 2, 3> b;
  b << 0.7, -0.2, 0.4,
      -0.3, 0.5, 0.1;
  Eigen::Matrix2d c;
  c << 2.1, 0.35,
       0.35, 1.4;
  Eigen::Matrix3d metric;
  metric << 1.6, 0.1, 0.0,
            0.1, 0.9, 0.12,
            0.0, 0.12, 1.2;
  const Eigen::Vector3d gradient(0.9, -1.1, 0.6);
  constexpr double radius = 0.35;

  NeoOptions options;
  options.trust_radius = radius;
  options.relative_residual_tolerance = 1.0e-12;
  const ResponseNeoResult result = xmvb::vb::solve_response_neo(
      dense_problem(a, b, c, metric, gradient), options);
  verify_coupled_residual(a, b, c, metric, gradient, radius, result);

  const Eigen::Matrix3d relaxed =
      a - b.transpose() * c.ldlt().solve(b);
  const Eigen::Vector3d reference = explicit_orbital_step(
      relaxed, metric, gradient, radius);
  require_close(result.step.orbital, reference, 2.0e-9,
                "projected-response step differs from explicit Schur NEO");
  require_close(result.step.structure, -c.ldlt().solve(b * reference),
                2.0e-9,
                "reconstructed structure response is inconsistent");
}

void check_anisotropic_coupled_dual_residual() {
  Eigen::Matrix2d a = Eigen::Matrix2d::Zero();
  a.diagonal() << 1.0, 1.0e4;
  const Eigen::RowVector2d b = Eigen::RowVector2d::Zero();
  const Eigen::Matrix<double, 1, 1> c =
      Eigen::Matrix<double, 1, 1>::Identity();
  Eigen::Matrix2d metric = Eigen::Matrix2d::Zero();
  metric.diagonal() << 1.0e-4, 1.0;
  const Eigen::Matrix2d metric_inverse =
      metric.ldlt().solve(Eigen::Matrix2d::Identity());
  const Eigen::Vector2d gradient(1.0e-4, 1.0e-4);
  const ResponseNeoProblem problem(
      gradient,
      1,
      [a, b](const Eigen::VectorXd& p) {
        return ResponseNeoDirection{a * p, b * p};
      },
      [b, c](const Eigen::VectorXd& q) {
        return ResponseNeoDirection{b.transpose() * q, c * q};
      },
      [metric](const Eigen::VectorXd& p) { return metric * p; },
      [metric_inverse](const Eigen::VectorXd& covector) {
        return std::sqrt(covector.dot(metric_inverse * covector));
      });

  NeoOptions options;
  options.trust_radius = 1.0;
  options.absolute_residual_tolerance = 1.0e-3;
  options.relative_residual_tolerance = 1.0e-8;
  options.require_curvature_certificate = false;
  const ResponseNeoResult result = xmvb::vb::solve_response_neo(
      problem, options);

  require(result.converged(), "anisotropic response NEO did not converge");
  require(result.orbital_actions >= 2,
          "anisotropic response NEO accepted its incomplete first direction");
  require_close(result.step.orbital,
                Eigen::Vector2d(-1.0e-4, -1.0e-8), 1.0e-10,
                "anisotropic response NEO missed the full Newton step");
  require(result.step.structure.stableNorm() <= 1.0e-12,
          "uncoupled anisotropic structure response was nonzero");
  const Eigen::Vector2d orbital_residual =
      gradient + a * result.step.orbital;
  const Eigen::VectorXd structure_residual =
      b * result.step.orbital + c * result.step.structure;
  const double physical_residual = std::hypot(
      physical_dual_norm(metric, orbital_residual),
      structure_residual.stableNorm());
  require(physical_residual <= result.residual_target,
          "anisotropic response NEO physical residual exceeds its target");
  require(std::abs(result.residual_norm - physical_residual) <= 1.0e-12,
          "anisotropic response NEO reported a coordinate residual norm");
}

void check_negative_relaxed_curvature() {
  Eigen::Matrix2d a;
  a << 0.4, 0.1,
       0.1, 1.8;
  Eigen::RowVector2d b;
  b << 1.2, 0.0;
  Eigen::Matrix<double, 1, 1> c;
  c << 1.0;
  const Eigen::Matrix2d metric = Eigen::Matrix2d::Identity();
  const Eigen::Vector2d gradient(0.3, -0.2);

  NeoOptions options;
  options.trust_radius = 0.5;
  options.relative_residual_tolerance = 1.0e-12;
  const ResponseNeoResult result = xmvb::vb::solve_response_neo(
      dense_problem(a, b, c, metric, gradient), options);
  verify_coupled_residual(
      a, b, c, metric, gradient, options.trust_radius, result);
  require(result.shift > 0.0,
          "negative relaxed curvature did not produce a spectral shift");
  require(result.boundary,
          "negative relaxed curvature did not reach the orbital boundary");
}

void check_budget_reports_subspace_limit() {
  const Eigen::Matrix3d a = Eigen::Matrix3d::Identity();
  Eigen::Matrix<double, 2, 3> b;
  b << 0.4, 0.2, 0.1,
      -0.1, 0.3, 0.5;
  Eigen::Matrix2d c;
  c << 1.5, 0.2,
       0.2, 1.0;
  const Eigen::Vector3d gradient(1.0, -0.4, 0.7);
  NeoOptions options;
  options.trust_radius = 0.2;
  options.relative_residual_tolerance = 1.0e-13;
  options.maximum_subspace_dimension = 1;
  const ResponseNeoResult result = xmvb::vb::solve_response_neo(
      dense_problem(
          a, b, c, Eigen::Matrix3d::Identity(), gradient), options);
  require(result.stop_reason == NeoStopReason::SubspaceLimit,
          "explicit response NEO budget was not reported");
}

void check_structure_contracts() {
  const Eigen::Matrix2d a = Eigen::Matrix2d::Identity();
  const Eigen::RowVector2d b(0.5, 0.0);
  const Eigen::Vector2d gradient(1.0, 0.2);
  NeoOptions options;
  options.trust_radius = 0.2;

  Eigen::Matrix<double, 1, 1> negative_c;
  negative_c << -1.0;
  const ResponseNeoResult excited_response = xmvb::vb::solve_response_neo(
      dense_problem(
          a, b, negative_c, Eigen::Matrix2d::Identity(), gradient),
      options);
  verify_coupled_residual(
      a, b, negative_c, Eigen::Matrix2d::Identity(), gradient,
      options.trust_radius, excited_response);

  Eigen::Matrix<double, 1, 1> singular_c;
  singular_c.setZero();
  bool rejected_range = false;
  try {
    (void)xmvb::vb::solve_response_neo(
        dense_problem(
            a, b, singular_c, Eigen::Matrix2d::Identity(), gradient),
        options);
  } catch (const std::runtime_error&) {
    rejected_range = true;
  }
  require(rejected_range, "coupling outside Range(C) was accepted");
}

void check_workspace_reuses_actions_after_radius_change() {
  Eigen::Matrix3d a;
  a << 2.0, 0.2, -0.1,
       0.2, 1.4, 0.3,
      -0.1, 0.3, 0.8;
  Eigen::Matrix<double, 2, 3> b;
  b << 0.4, -0.3, 0.2,
       0.1, 0.5, -0.2;
  Eigen::Matrix2d c;
  c << 1.8, 0.2,
       0.2, 1.1;
  const Eigen::Matrix3d metric = Eigen::Matrix3d::Identity();
  const Eigen::Vector3d gradient(0.7, -0.8, 0.5);
  int gradient_norm_calls = 0;
  const ResponseNeoProblem problem(
      gradient, c.rows(),
      [a, b](const Eigen::VectorXd& p) {
        return ResponseNeoDirection{a * p, b * p};
      },
      [b, c](const Eigen::VectorXd& q) {
        return ResponseNeoDirection{b.transpose() * q, c * q};
      },
      [metric](const Eigen::VectorXd& p) { return metric * p; },
      [&](const Eigen::VectorXd& covector) {
        if ((covector - gradient).isZero(0.0)) ++gradient_norm_calls;
        return physical_dual_norm(metric, covector);
      });
  ResponseNeoWorkspace workspace(problem);

  NeoOptions options;
  options.trust_radius = 0.5;
  options.relative_residual_tolerance = 1.0e-12;
  const ResponseNeoResult first = workspace.solve(options);
  require(first.converged(), "initial reusable NEO solve did not converge");
  require(first.coupled_actions > 0,
          "initial reusable NEO solve performed no actions");

  options.trust_radius = 0.2;
  const ResponseNeoResult reused = workspace.solve(options);
  require(gradient_norm_calls == 1,
          "radius retry recomputed the accepted-point gradient dual norm");
  const ResponseNeoResult fresh = xmvb::vb::solve_response_neo(
      problem, options);
  require(reused.converged(), "reused NEO solve did not converge");
  require(fresh.converged(), "fresh comparison NEO solve did not converge");
  require(reused.coupled_actions == 0,
          "radius-only retry repeated cached coupled actions");
  require(reused.coupled_actions < fresh.coupled_actions,
          "reused NEO solve did not reduce coupled actions");
  require_close(reused.step.orbital, fresh.step.orbital, 2.0e-10,
                "reused NEO orbital step differs from a fresh solve");
  require_close(reused.step.structure, fresh.step.structure, 2.0e-10,
                "reused NEO response differs from a fresh solve");
}

void check_recycled_orbital_guess_starts_subspace() {
  Eigen::VectorXd first_direction;
  const Eigen::Vector2d guess(0.0, 2.0);
  const ResponseNeoProblem problem(
      Eigen::Vector2d(1.0, 0.3),
      0,
      [&first_direction](const Eigen::VectorXd& p) {
        if (first_direction.size() == 0) first_direction = p;
        return ResponseNeoDirection{p, Eigen::VectorXd::Zero(0)};
      },
      [](const Eigen::VectorXd&) {
        return ResponseNeoDirection{
            Eigen::VectorXd::Zero(2), Eigen::VectorXd::Zero(0)};
      },
      [](const Eigen::VectorXd& p) { return p; },
      [](const Eigen::VectorXd& covector) { return covector.stableNorm(); },
      {},
      guess);
  NeoOptions options;
  options.trust_radius = 0.5;
  options.relative_residual_tolerance = 1.0e-12;
  const ResponseNeoResult result = xmvb::vb::solve_response_neo(
      problem, options);
  require(result.converged(), "recycled-guess NEO solve did not converge");
  require_close(first_direction, guess.normalized(), 1.0e-14,
                "recycled orbital guess did not start the Davidson space");
}

void check_boundary_residual_uses_shifted_preconditioner() {
  double largest_shift = 0.0;
  Eigen::Matrix3d hessian;
  hessian << -1.2, 0.2, 0.1,
              0.2, 0.8, 0.3,
              0.1, 0.3, 1.7;
  const ResponseNeoProblem problem(
      Eigen::Vector3d(0.4, -0.7, 0.2),
      0,
      [hessian](const Eigen::VectorXd& p) {
        return ResponseNeoDirection{
            hessian * p, Eigen::VectorXd::Zero(0)};
      },
      [](const Eigen::VectorXd&) {
        return ResponseNeoDirection{
            Eigen::VectorXd::Zero(3), Eigen::VectorXd::Zero(0)};
      },
      [](const Eigen::VectorXd& p) { return p; },
      [](const Eigen::VectorXd& covector) { return covector.stableNorm(); },
      [&largest_shift](
          const Eigen::VectorXd& residual,
          double shift,
          double) {
        largest_shift = std::max(largest_shift, shift);
        return residual / (2.0 + shift);
      });
  NeoOptions options;
  options.trust_radius = 0.2;
  options.relative_residual_tolerance = 1.0e-12;
  const ResponseNeoResult result = xmvb::vb::solve_response_neo(
      problem, options);
  require(result.converged(), "shift-aware response NEO did not converge");
  require(result.boundary && result.shift > 0.0,
          "shift-aware response NEO fixture did not reach the boundary");
  require(largest_shift > 0.0,
          "boundary KKT residual did not pass its shift to the preconditioner");
}

void check_preconditioner_receives_relative_physical_budget() {
  const Eigen::Vector2d gradient(1.0, 1.0);
  const Eigen::Vector2d metric_diagonal(1.0e-4, 1.0e4);
  const double gradient_norm = std::sqrt(1.0e4 + 1.0e-4);
  double observed_fraction = -1.0;
  const ResponseNeoProblem problem(
      gradient, 0,
      [](const Eigen::VectorXd& p) {
        return ResponseNeoDirection{p, Eigen::VectorXd::Zero(0)};
      },
      [](const Eigen::VectorXd&) {
        return ResponseNeoDirection{
            Eigen::VectorXd::Zero(2), Eigen::VectorXd::Zero(0)};
      },
      [metric_diagonal](const Eigen::VectorXd& p) {
        return (metric_diagonal.array() * p.array()).matrix().eval();
      },
      [metric_diagonal](const Eigen::VectorXd& covector) {
        return std::sqrt((covector.array().square() /
                          metric_diagonal.array()).sum());
      },
      [&observed_fraction](const Eigen::VectorXd& residual,
                           double, double fraction) {
        observed_fraction = fraction;
        return residual;
      });
  NeoOptions options;
  options.trust_radius = 1.0;
  options.relative_residual_tolerance = 0.25;
  options.absolute_residual_tolerance = 0.1;
  options.maximum_subspace_dimension = 1;
  options.require_curvature_certificate = false;
  static_cast<void>(xmvb::vb::solve_response_neo(problem, options));
  const double expected_fraction =
      0.25 + options.absolute_residual_tolerance / gradient_norm;
  require(std::abs(observed_fraction - expected_fraction) <= 1.0e-12,
          "orbital preconditioner received an absolute physical target");
}

void check_inexact_block_preconditioner_preserves_coupled_step() {
  Eigen::Matrix3d a;
  a << 3.0, 0.3, -0.1,
       0.3, 2.1, 0.2,
       -0.1, 0.2, 1.4;
  Eigen::RowVector3d b;
  b << 0.6, -0.2, 0.3;
  Eigen::Matrix<double, 1, 1> c;
  c << 1.7;
  Eigen::Matrix3d metric;
  metric << 1.5, 0.2, 0.0,
            0.2, 1.1, -0.1,
            0.0, -0.1, 0.8;
  const Eigen::Vector3d gradient(0.9, -0.7, 0.5);
  const Eigen::Vector3d block_diagonal(3.0, 2.1, 1.4);
  const Eigen::Matrix3d metric_inverse =
      metric.ldlt().solve(Eigen::Matrix3d::Identity());
  double largest_shift = 0.0;
  const ResponseNeoProblem problem(
      gradient,
      1,
      [a, b](const Eigen::VectorXd& p) {
        return ResponseNeoDirection{a * p, b * p};
      },
      [b, c](const Eigen::VectorXd& q) {
        return ResponseNeoDirection{b.transpose() * q, c * q};
      },
      [metric](const Eigen::VectorXd& p) { return metric * p; },
      [metric_inverse](const Eigen::VectorXd& covector) {
        return std::sqrt(covector.dot(metric_inverse * covector));
      },
      [&largest_shift, block_diagonal](
          const Eigen::VectorXd& residual,
          double shift,
          double residual_target) {
        largest_shift = std::max(largest_shift, shift);
        require(residual_target > 0.0,
                "NEO did not pass its KKT target to the preconditioner");
        // Deliberately omit the off-diagonal physical metric: this is only
        // a cheap positive preconditioner, not the shifted Newton inverse.
        return (residual.array() /
                (block_diagonal.array() + shift)).matrix();
      });
  NeoOptions options;
  options.trust_radius = 0.25;
  options.relative_residual_tolerance = 1.0e-12;
  const ResponseNeoResult result = xmvb::vb::solve_response_neo(
      problem, options);
  verify_coupled_residual(
      a, b, c, metric, gradient, options.trust_radius, result);
  const Eigen::Matrix3d relaxed =
      a - b.transpose() * c.ldlt().solve(b);
  require_close(
      result.step.orbital,
      explicit_orbital_step(
          relaxed, metric, gradient, options.trust_radius),
      2.0e-9,
      "inexact block preconditioner changed the coupled NEO step");
  require(largest_shift > 0.0,
          "inexact block preconditioner did not receive the NEO shift");
}

void check_keyframe_response_uses_forcing_accuracy() {
  Eigen::Matrix4d a = Eigen::Matrix4d::Zero();
  a.diagonal() << 2.0, 3.0, 4.0, 5.0;
  Eigen::Matrix<double, 3, 4> b =
      Eigen::Matrix<double, 3, 4>::Zero();
  b.diagonal() << 0.5, 0.4, 0.3;
  Eigen::Matrix3d c = Eigen::Matrix3d::Zero();
  c.diagonal() << 2.0, 1.5, 1.0;
  const Eigen::Vector4d gradient(1.0, 0.0, 0.0, 0.0);
  const ResponseNeoProblem problem = dense_problem(
      a, b, c, Eigen::Matrix4d::Identity(), gradient);

  NeoOptions keyframe_options;
  keyframe_options.trust_radius = 1.0;
  keyframe_options.relative_residual_tolerance = 0.1;
  keyframe_options.absolute_residual_tolerance = 1.0e-3;
  keyframe_options.require_curvature_certificate = false;
  const ResponseNeoResult keyframe = xmvb::vb::solve_response_neo(
      problem, keyframe_options);
  require(keyframe.converged(),
          "inexact keyframe response solve did not converge");
  require(keyframe.minimum_curvature_orbital.norm() > 0.0 &&
              std::isfinite(keyframe.curvature_residual_norm) &&
              keyframe.curvature_residual_target > 0.0,
          "deferred curvature diagnostic was not evaluated at convergence");
  require(problem.orbital_dual_norm(keyframe.kkt_residual.orbital) <=
              keyframe.residual_target,
          "inexact keyframe orbital residual exceeds its certificate");
  require(keyframe.kkt_residual.structure.stableNorm() <=
              keyframe.residual_target,
          "inexact keyframe structure residual exceeds its certificate");

  NeoOptions certified_options = keyframe_options;
  certified_options.require_curvature_certificate = true;
  const ResponseNeoResult certified = xmvb::vb::solve_response_neo(
      problem, certified_options);
  require(certified.converged(),
          "curvature-certified comparison solve did not converge");
  require(keyframe.structure_actions <= certified.structure_actions,
          "keyframe forcing accuracy increased structure actions");
}

void check_response_is_local_to_candidate_direction() {
  constexpr int size = 5;
  Eigen::MatrixXd a = Eigen::MatrixXd::Zero(size, size);
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(size, size);
  Eigen::MatrixXd c = Eigen::MatrixXd::Zero(size, size);
  a.diagonal().setLinSpaced(2.0, 6.0);
  b.diagonal().setLinSpaced(0.2, 0.6);
  c.diagonal().setLinSpaced(1.0, 2.0);
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(size);
  gradient[0] = 1.0;

  NeoOptions options;
  options.trust_radius = 1.0;
  options.relative_residual_tolerance = 1.0e-12;
  const ResponseNeoResult result = xmvb::vb::solve_response_neo(
      dense_problem(a, b, c, Eigen::MatrixXd::Identity(size, size), gradient),
      options);
  verify_coupled_residual(
      a, b, c, Eigen::MatrixXd::Identity(size, size), gradient,
      options.trust_radius, result);
  require(result.structure_actions == 1,
          "response NEO solved unused orbital-basis response directions");
}

void check_small_gradient_needs_relative_accuracy() {
  Eigen::Matrix3d a = Eigen::Vector3d(1.0, 3.0, 9.0).asDiagonal();
  Eigen::Matrix3d b = 0.2 * Eigen::Matrix3d::Identity();
  Eigen::Matrix3d c = 2.0 * Eigen::Matrix3d::Identity();
  const Eigen::Vector3d gradient(1.0e-5, -2.0e-5, 3.0e-5);
  const auto problem = dense_problem(a, b, c, Eigen::Matrix3d::Identity(), gradient);
  NeoOptions options;
  options.trust_radius = 1.0;
  options.relative_residual_tolerance = gradient.norm() / (1.0 + gradient.norm());
  options.require_curvature_certificate = false;
  const auto result = xmvb::vb::solve_response_neo(problem, options);
  const Eigen::Matrix3d relaxed = a - b.transpose() * c.ldlt().solve(b);
  const Eigen::Vector3d reference = -relaxed.ldlt().solve(gradient);
  require(result.converged(), "small-gradient relative solve did not converge");
  const double target = options.relative_residual_tolerance *
      problem.orbital_dual_norm(gradient);
  // Eliminating the structure residual gives H_relaxed (p-p_exact)
  // = r_o - B^T C^-1 r_s. Bound the step error by the requested KKT accuracy,
  // rather than demanding accuracy stricter than the solver was asked for.
  const double error_bound = relaxed.ldlt().solve(Eigen::Matrix3d::Identity()).norm() *
      std::hypot(1.0, c.ldlt().solve(b).norm()) * target;
  require((result.step.orbital - reference).norm() <= error_bound +
              64.0 * std::numeric_limits<double>::epsilon() * reference.norm(),
          "small-gradient response step exceeds its KKT-derived error bound");
  require(result.residual_norm <= target,
          "small-gradient KKT residual was hidden by an absolute tolerance");
}

void check_small_structure_residual_requires_response_correction() {
  Eigen::Matrix2d a = Eigen::Matrix2d::Zero();
  a.diagonal() << 2.0, 3.0;
  const Eigen::RowVector2d b(1.0e-3, 0.0);
  const Eigen::Matrix<double, 1, 1> c =
      1.0e-6 * Eigen::Matrix<double, 1, 1>::Identity();
  const Eigen::Matrix2d metric = Eigen::Matrix2d::Identity();
  const Eigen::Vector2d gradient(0.1, 0.0);
  int callback_calls = 0;
  int structure_image_calls = 0;
  double first_structure_residual = 0.0;
  const ResponseNeoProblem problem(
      gradient,
      1,
      [a, b](const Eigen::VectorXd& p) {
        return ResponseNeoDirection{a * p, b * p};
      },
      [b, c, &structure_image_calls](const Eigen::VectorXd& q) {
        ++structure_image_calls;
        return ResponseNeoDirection{b.transpose() * q, c * q};
      },
      [metric](const Eigen::VectorXd& p) { return metric * p; },
      [](const Eigen::VectorXd& covector) { return covector.stableNorm(); },
      {}, {}, 0.0, {},
      [&callback_calls, &first_structure_residual, c](
          const Eigen::VectorXd& residual) {
        if (callback_calls == 0) {
          first_structure_residual = residual.stableNorm();
        }
        ++callback_calls;
        return (-residual / c(0, 0)).eval();
      });
  ResponseNeoWorkspace workspace(problem);
  NeoOptions options;
  options.trust_radius = 1.0;
  options.relative_residual_tolerance = 1.0e-12;
  options.absolute_residual_tolerance = 1.0e-3;
  options.require_curvature_certificate = false;
  const ResponseNeoResult result = workspace.solve(options);

  require(result.converged(), "response-corrected NEO did not converge");
  require(callback_calls >= 1 &&
              first_structure_residual > 0.0 &&
              first_structure_residual < options.absolute_residual_tolerance,
          "small structure residual did not trigger a response solve");
  require(workspace.structure_basis_size() == 1,
          "missing response direction was not retained");
  require(structure_image_calls == callback_calls,
          "retaining the response correction repeated its coupled action: " +
              std::to_string(structure_image_calls) + " images for " +
              std::to_string(callback_calls) + " residual solves");
  const Eigen::Matrix2d relaxed =
      a - b.transpose() * (b / c(0, 0));
  require_close(result.step.orbital,
                explicit_orbital_step(relaxed, metric, gradient,
                                      options.trust_radius),
                2.0e-9,
                "small structure residual hid a large orbital response");
  verify_coupled_residual(a, b, c, metric, gradient,
                          options.trust_radius, result);
}

void check_incomplete_structure_basis_cannot_certify_global_curvature() {
  Eigen::Matrix2d a = Eigen::Matrix2d::Zero();
  a.diagonal() << -1.0, 1.0;
  const Eigen::RowVector2d b(0.0, 3.0);
  const Eigen::Matrix<double, 1, 1> c =
      Eigen::Matrix<double, 1, 1>::Identity();
  const Eigen::Matrix2d metric = Eigen::Matrix2d::Identity();
  const Eigen::Vector2d gradient(0.1, 0.0);
  const ResponseNeoProblem problem = dense_problem(a, b, c, metric, gradient);
  ResponseNeoWorkspace workspace(problem);

  NeoOptions options;
  options.trust_radius = 0.5;
  options.relative_residual_tolerance = 1.0e-12;
  options.require_curvature_certificate = true;
  const ResponseNeoResult result = workspace.solve(options);
  require(result.converged(), "partial-response fixture did not converge");
  require(workspace.orbital_basis_size() == 2 &&
              workspace.structure_basis_size() == 0,
          "partial-response fixture did not retain only the orbital basis");
  require(result.curvature_residual_norm <= result.curvature_residual_target,
          "partial-response lowest Ritz pair was not converged");
  // The retained lowest mode is -1, but the omitted response to the second
  // orbital lowers its relaxed curvature from +1 to 1 - 3^2 = -8.
  const Eigen::Matrix2d relaxed =
      a - b.transpose() * (b / c(0, 0));
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> full_spectrum(relaxed);
  require(full_spectrum.info() == Eigen::Success,
          "partial-response full spectrum failed");
  const double retained_curvature =
      result.minimum_curvature_orbital.dot(
          a * result.minimum_curvature_orbital) /
      result.minimum_curvature_orbital.squaredNorm();
  require(full_spectrum.eigenvalues()[0] < retained_curvature - 1.0,
          "partial-response fixture has no hidden lower curvature");
  require(!result.global_curvature_certified,
          "incomplete structure response falsely certified global curvature");
}

void check_early_trial_is_not_certified_and_resumes() {
  Eigen::Matrix<double, 6, 6> a = Eigen::Matrix<double, 6, 6>::Zero();
  a.diagonal() << 1.0, 2.0, 4.0, 8.0, 16.0, 32.0;
  const Eigen::MatrixXd b(0, 6), c(0, 0);
  const auto metric = Eigen::Matrix<double, 6, 6>::Identity().eval();
  const Eigen::VectorXd gradient = 0.01 * Eigen::VectorXd::Ones(6);
  const auto problem = dense_problem(a, b, c, metric, gradient);
  ResponseNeoWorkspace workspace(problem);
  NeoOptions options;
  options.trust_radius = 0.001;
  options.relative_residual_tolerance = 1.0e-12;
  options.require_curvature_certificate = false;
  const auto early = workspace.solve(options, true);
  require(early.trial_ready() && !early.converged(),
          "early trial falsely reported Newton convergence");
  require(early.residual_norm > early.residual_target,
          "early-trial fixture already satisfies the strict certificate");
  require(workspace.orbital_basis_size() == 3,
          "CIAH trial did not start after three orbital directions");
  int retained = static_cast<int>(workspace.orbital_basis_size());
  options.trust_radius *= 0.5;
  const auto retry = workspace.solve(options, true);
  require(workspace.orbital_basis_size() > retained,
          "rejected early trial yielded again without refining its subspace");
  require(retry.orbital_actions == workspace.orbital_basis_size() - retained,
          "early radius retry repeated retained HVPs");
  retained = static_cast<int>(workspace.orbital_basis_size());
  const auto resumed = workspace.solve(options);
  verify_coupled_residual(a, b, c, metric, gradient, options.trust_radius, resumed);
  require(resumed.orbital_actions == workspace.orbital_basis_size() - retained,
          "resuming an early trial repeated retained HVPs");
  require_close(resumed.step.orbital,
      explicit_orbital_step(a, metric, gradient, options.trust_radius),
      2.0e-9, "resumed early trial changed the strict trust-region solution");
}

void check_historical_interior_trial_preserves_physical_certificate() {
  const Eigen::MatrixXd a = Eigen::Vector4d(1.0, 2.0, 4.0, 8.0).asDiagonal();
  const Eigen::MatrixXd b(0, 4), c(0, 0);
  const Eigen::MatrixXd metric = 1.0e-4 * Eigen::MatrixXd::Identity(4, 4);
  const Eigen::VectorXd gradient = Eigen::VectorXd::Ones(4);
  const auto problem = dense_problem(a, b, c, metric, gradient);
  ResponseNeoWorkspace workspace(problem);
  NeoOptions options;
  options.trust_radius = 10.0;
  options.relative_residual_tolerance = 1.0e-12;
  options.require_curvature_certificate = false;
  const auto trial = workspace.solve(options, true);
  require(trial.trial_ready() && !trial.converged() && !trial.boundary,
          "historical interior candidate was not classified as an uncertified trial");
  require(trial.residual_norm > trial.residual_target,
          "historical trial fixture already satisfies the physical certificate");
  const int retained = static_cast<int>(workspace.orbital_basis_size());
  const auto final = workspace.solve(options);
  verify_coupled_residual(a, b, c, metric, gradient, options.trust_radius, final);
  require(final.orbital_actions == workspace.orbital_basis_size() - retained,
          "strict continuation recomputed historical trial HVPs");
  require_close(final.step.orbital, -a.ldlt().solve(gradient), 1.0e-10,
                "historical trial changed the certified Newton solution");
}

void check_trial_defers_structure_certificate() {
  int response_calls = 0;
  const ResponseNeoProblem problem(
      Eigen::VectorXd::Ones(1), 1,
      [](const Eigen::VectorXd& p) {
        return ResponseNeoDirection{2 * p, 0.1 * p};
      },
      [](const Eigen::VectorXd& q) {
        return ResponseNeoDirection{0.1 * q, q};
      },
      [](const Eigen::VectorXd& p) { return p; },
      [](const Eigen::VectorXd& r) { return r.norm(); },
      {}, {}, 0.0, {},
      [&response_calls](const Eigen::VectorXd& r) {
        ++response_calls;
        return Eigen::VectorXd(-r);
      });
  NeoOptions options;
  options.relative_residual_tolerance = 0.1;
  options.require_curvature_certificate = false;
  ResponseNeoWorkspace workspace(problem);
  const auto trial = workspace.solve(options, true);
  require(trial.trial_ready() && !trial.converged() && response_calls == 0,
          "nonlinear trial paid for or claimed exact structure certification");
  const auto certified = workspace.solve(options);
  require(certified.converged() && response_calls == 1 &&
              certified.residual_norm <= certified.residual_target,
          "strict continuation bypassed the structure response certificate");
  require_close(certified.kkt_residual.structure,
                0.1 * certified.step.orbital + certified.step.structure, 1e-12,
                "deferred certificate returned an inconsistent structure residual");
}

void check_physical_mismatch_seeds_response() {
  Eigen::Matrix2d a;
  a << 2.0, 0.6, 0.6, 3.0;
  const Eigen::RowVector2d b(0.2, 0.1);
  const Eigen::Matrix<double, 1, 1> c =
      1.7 * Eigen::Matrix<double, 1, 1>::Identity();
  const Eigen::Matrix2d metric = 1.0e-4 * Eigen::Matrix2d::Identity();
  const Eigen::Vector2d gradient(0.0005, 0.001);
  int orbital_actions = 0, responses = 0;
  const ResponseNeoProblem problem(
      gradient, 1,
      [&](const Eigen::VectorXd& p) {
        ++orbital_actions;
        return ResponseNeoDirection{a * p, b * p};
      },
      [b, c](const Eigen::VectorXd& q) {
        return ResponseNeoDirection{b.transpose() * q, c * q};
      },
      [metric](const Eigen::VectorXd& p) { return metric * p; },
      [](const Eigen::VectorXd& r) { return 100.0 * r.norm(); },
      {}, {}, 0.0, {},
      [&](const Eigen::VectorXd& r) {
        if (responses == 0)
          require(orbital_actions == 1,
                  "physical mismatch accumulated frozen orbital directions");
        ++responses;
        return Eigen::VectorXd(-r / c(0, 0));
      });
  NeoOptions options;
  options.trust_radius = 1.0;
  options.relative_residual_tolerance = 1.0e-12;
  options.absolute_residual_tolerance = 1.0e-3;
  options.require_curvature_certificate = false;
  ResponseNeoWorkspace workspace(problem);
  const auto result = workspace.solve(options);
  require(responses > 0, "physical mismatch omitted its structure seed");
  verify_coupled_residual(a, b, c, metric, gradient, options.trust_radius, result);
  const Eigen::Matrix2d relaxed = a - b.transpose() * c.ldlt().solve(b);
  require_close(result.step.orbital, -relaxed.ldlt().solve(gradient), 1.0e-10,
                "direction-local seeding changed the coupled Newton solution");
}

void check_candidate_observation_preserves_solve() {
  const Eigen::MatrixXd a = Eigen::Vector4d(1.0, 2.0, 4.0, 8.0).asDiagonal();
  const Eigen::MatrixXd b(0, 4), c(0, 0);
  const Eigen::MatrixXd metric = Eigen::MatrixXd::Identity(4, 4);
  const Eigen::VectorXd gradient = Eigen::VectorXd::Ones(4);
  const auto problem = dense_problem(a, b, c, metric, gradient);
  NeoOptions options;
  options.trust_radius = 10.0;
  options.relative_residual_tolerance = 1.0e-12;
  options.absolute_residual_tolerance = 1.0e-12;
  options.require_curvature_certificate = false;
  ResponseNeoWorkspace reference_workspace(problem), observed_workspace(problem);
  const auto reference = reference_workspace.solve(options);
  int observations = 0;
  int model_observations = 0;
  const auto observed = observed_workspace.solve(options, false,
      [&](const ResponseNeoResult& candidate) {
        ++observations;
        require(candidate.step.orbital.allFinite(), "observed step is not finite");
        require(candidate.predicted_reduction > 0.0, "observed model has no descent");
        require(!candidate.converged(), "pre-certificate observation claimed convergence");
      },
      [&](const Eigen::MatrixXd& basis, const Eigen::MatrixXd& core,
          const Eigen::MatrixXd& model) {
        ++model_observations;
        require_close(basis.transpose() * basis,
                      Eigen::MatrixXd::Identity(basis.cols(), basis.cols()),
                      1.0e-13, "observed model basis is not orthonormal");
        require_close(model, basis.transpose() * a * basis, 1.0e-13,
                      "observed projected model disagrees with the dense action");
        require_close(core, model, 1.0e-13,
                      "zero structure space changed the observed relaxed model");
      });
  require(observations == observed.iterations, "not every projection was observed");
  require(model_observations == observations, "not every new model was observed");
  require(observed.orbital_actions == reference.orbital_actions,
          "observation added HVP actions");
  require_close(observed.step.orbital, reference.step.orbital, 1.0e-13,
                "observation changed the final step");
  verify_coupled_residual(a, b, c, metric, gradient, options.trust_radius, observed);
}

}  // namespace

int main() {
  try {
    check_physical_mismatch_seeds_response();
    check_candidate_observation_preserves_solve();
    check_early_trial_is_not_certified_and_resumes();
    check_historical_interior_trial_preserves_physical_certificate();
    check_trial_defers_structure_certificate();
    check_spd_structure_and_explicit_schur();
    check_anisotropic_coupled_dual_residual();
    check_negative_relaxed_curvature();
    check_budget_reports_subspace_limit();
    check_structure_contracts();
    check_workspace_reuses_actions_after_radius_change();
    check_recycled_orbital_guess_starts_subspace();
    check_boundary_residual_uses_shifted_preconditioner();
    check_preconditioner_receives_relative_physical_budget();
    check_inexact_block_preconditioner_preserves_coupled_step();
    check_keyframe_response_uses_forcing_accuracy();
    check_response_is_local_to_candidate_direction();
    check_small_gradient_needs_relative_accuracy();
    check_small_structure_residual_requires_response_correction();
    check_incomplete_structure_basis_cannot_certify_global_curvature();
    std::cout << "response NEO solver tests passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "response NEO solver test failed: " << error.what() << '\n';
    return 1;
  }
}
