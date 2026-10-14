#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include "vbscf/optimization/neo/ciah.hpp"
#include "vbscf/optimization/backends/ciah.hpp"
#include "vbscf/optimization/preconditioners/coarse.hpp"

namespace {

using xmvb::vb::CiahOptions;
using xmvb::vb::CiahOriginRay;
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
    const Eigen::VectorXd& gradient, int* actions = nullptr,
    int* metrics = nullptr, int* norms = nullptr) {
  return NeoProblem(
      gradient,
      [hessian, actions](const Eigen::VectorXd& vector) {
        if (actions) ++*actions;
        return Eigen::VectorXd(hessian * vector);
      },
      [metric, metrics](const Eigen::VectorXd& vector) {
        if (metrics) ++*metrics;
        return Eigen::VectorXd(metric * vector);
      },
      [metric, norms](const Eigen::VectorXd& covector) {
        if (norms) ++*norms;
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
  const Eigen::VectorXd expected = result.reference_component *
      (gradient + hessian * result.step -
       result.eigenvalue * metric * result.step);
  // The eigenvector phase is free, so the stored residual may have either sign.
  const double error = std::min(
      (result.residual - expected).stableNorm(),
      (result.residual + expected).stableNorm());
  require(error <= 2.0e-12 * std::max(1.0, expected.stableNorm()),
          "CIAH returned an inconsistent full-space AH residual");
  require(std::abs(result.residual_norm -
                   result.residual.norm()) <= 2.0e-12,
          "CIAH residual norm omitted full coordinate components");
  const double expected_norm = (gradient + hessian * result.step).norm();
  require(std::abs(result.model_gradient_norm - expected_norm) <=
              2.0e-11 * std::max(1.0, expected_norm),
          "CIAH model gradient norm differs from the full coordinate norm");
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

/** @brief Recycled CI directions transform with a relaxed root's arbitrary phase. */
void check_restart_phase() {
  using xmvb::vb::optimizer_detail::CiahMacroStep;
  const Eigen::Vector3d reference(0.8,0.4,0.1);
  Eigen::Matrix3d metric;
  metric << 1.0,0.2,0.1,0.2,2.0,0.3,0.1,0.3,3.0;
  const Eigen::Vector3d root(0.7,0.5,0.2);
  const Eigen::Vector3d response(0.1,-0.2,0.05);
  CiahMacroStep positive,negative;
  positive.next_structure_guess=response;
  negative.next_structure_guess=response;
  positive.next_orbital_guess_packed=reference;
  negative.next_orbital_guess_packed=reference;
  positive.align_structure_phase(reference,metric*root);
  negative.align_structure_phase(reference,-metric*root);
  require_close(positive.next_structure_guess,response,1.0e-14,
      "CIAH changed a phase-consistent structure seed");
  require_close(negative.next_structure_guess,-response,1.0e-14,
      "CIAH failed to reverse the seed with the relaxed root");
  require_close(positive.next_orbital_guess_packed,negative.next_orbital_guess_packed,
      1.0e-14,"CIAH phase alignment changed the orbital direction");
  // B changes sign with the CI root; the coupled orbital image is invariant.
  Eigen::Matrix3d coupling;
  coupling << 2.0,0.3,0.1,0.2,1.0,0.4,0.3,-0.1,0.7;
  require_close(coupling*positive.next_structure_guess,
      -coupling*negative.next_structure_guess,1.0e-14,
      "CIAH restart violates coupled-action phase covariance");
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
  CiahOptions options;
  int observed = 0;
  options.observe_images = [&](const Eigen::VectorXd& v,const Eigen::VectorXd& gv) {
    ++observed;
    require_close(gv,metric*v,2.0e-12,"CIAH observer did not reuse the metric image");
    require(actions == observed,"CIAH observer is not attached to the matching HVP");
  };
  CiahWorkspace workspace(problem, Eigen::VectorXd::Ones(3),options);
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
  int joint_calls = 0;
  CiahOptions joint_options;
  joint_options.apply_images = [&](const Eigen::VectorXd& v) {
    ++joint_calls;
    return std::make_pair((hessian*v).eval(), (metric*v).eval());
  };
  const int separate_calls = actions;
  CiahWorkspace joint(problem, Eigen::VectorXd::Ones(3), joint_options);
  CiahResult joint_result;
  for (int i = 0; i < 5; ++i) {
    joint_result = joint.next(gradient);
    check_images_and_residual(hessian, metric, gradient, joint_result);
    if (joint_result.converged) break;
  }
  require(joint_result.converged && actions == separate_calls &&
      joint_calls == joint_result.basis_size,
      "joint CIAH action repeated a physical operator");
  require_close(joint_result.step, result.step, 2e-12,
      "joint H/G action changed the CIAH solution");
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

void check_origin_ray_residual(const CiahOriginRay& ray,
    const Eigen::MatrixXd& hessian, const Eigen::MatrixXd& metric,
    const Eigen::VectorXd& gradient) {
  require_close(ray.hessian_step, hessian * ray.step, 2.0e-12,
      "origin ray Hessian image does not match the cached action");
  const Eigen::VectorXd image = hessian * ray.step -
      ray.eigenvalue * metric * ray.step;
  const double expected = std::min(
      (image + ray.reference_component * gradient).stableNorm(),
      (image - ray.reference_component * gradient).stableNorm());
  require(std::abs(ray.residual_norm - expected) <
      2.0e-11 * std::max(1.0, expected),
      "origin ray full augmented residual is inconsistent");
}

void check_cached_origin_ray() {
  const Eigen::MatrixXd metric = Eigen::MatrixXd::Identity(1, 1);
  for (const double curvature : {2.0, -100.0}) {
    Eigen::MatrixXd hessian = Eigen::MatrixXd::Constant(1, 1, curvature);
    const Eigen::VectorXd gradient = Eigen::VectorXd::Constant(
        1, curvature > 0.0 ? 1.0 : 0.001);
    int h_calls = 0, g_calls = 0, full_norm_calls = 0;
    int preconditioner_calls = 0, coordinate_norm_calls = 0;
    const NeoProblem problem = dense_problem(
        hessian, metric, gradient, &h_calls, &g_calls, &full_norm_calls);
    CiahOptions options;
    options.coordinate_norm = [&](const Eigen::VectorXd& v) {
      ++coordinate_norm_calls;
      return v.norm();
    };
    options.residual_preconditioner =
        [&](const Eigen::VectorXd& v, double) {
          ++preconditioner_calls;
          return v;
        };
    CiahWorkspace workspace(problem, Eigen::VectorXd::Ones(1), options);
    const CiahResult usual = workspace.next(gradient);
    const int h_before = h_calls, g_before = g_calls;
    const int norm_before = full_norm_calls;
    const int preconditioner_before = preconditioner_calls;
    const int coordinate_before = coordinate_norm_calls;
    const auto ray = workspace.origin_ray(gradient);
    require(ray.has_value(), "cached origin AH has no ray in a one-dimensional span");
    require(h_calls == h_before && g_calls == g_before &&
        full_norm_calls == norm_before &&
        preconditioner_calls == preconditioner_before &&
        coordinate_norm_calls == coordinate_before + 1,
        "origin ray reevaluated an operator or physical norm, or missed its coordinate norm");
    require(gradient.dot(ray->step) < 0.0,
        "origin AH ray did not orient against its gradient");
    require(std::abs(ray->eigenvalue -
        0.5 * (curvature - std::sqrt(curvature * curvature +
                                     4.0 * gradient.squaredNorm()))) < 1.0e-11,
        "origin ray did not select the lowest AH eigenvalue");
    check_origin_ray_residual(*ray, hessian, metric, gradient);
    if (curvature < 0.0) {
      require(usual.eigenvalue > 0.0 && ray->eigenvalue < -99.0 &&
          ray->reference_component < 0.1,
          "weakly coupled negative mode did not bypass the usual root filter");
    }
  }

  Eigen::Matrix2d hessian = Eigen::Matrix2d::Zero();
  hessian.diagonal() << -100.0, 2.0;
  const Eigen::Matrix2d two_metric = Eigen::Matrix2d::Identity();
  const Eigen::Vector2d gradient(0.0, 1.0);
  int h_calls = 0, g_calls = 0, full_norm_calls = 0;
  const NeoProblem problem = dense_problem(
      hessian, two_metric, gradient, &h_calls, &g_calls, &full_norm_calls);
  CiahWorkspace workspace(problem, Eigen::VectorXd::Unit(2, 0));
  const CiahResult usual = workspace.next(gradient);
  require(usual.basis_size == 1 && usual.eigenvalue == 0.0,
      "zero-reference fixture unexpectedly captured the full gradient");
  const auto ray = workspace.origin_ray(gradient);
  require(ray && ray->reference_component < 1.0e-14 &&
      gradient.dot(ray->step) == 0.0 &&
      ray->step.dot(ray->hessian_step) < 0.0,
      "zero-reference negative-curvature tail was discarded");
  check_origin_ray_residual(*ray, hessian, two_metric, gradient);
  require(workspace.projected_dual_norm(gradient) == 0.0 &&
      gradient.norm() > 0.0 && h_calls == 1 && g_calls == 1 &&
      full_norm_calls == 0,
      "zero-reference ray falsely captured the gradient or repeated an action");

  const Eigen::Matrix2d positive_hessian =
      Eigen::Vector2d(2.0, 3.0).asDiagonal();
  const Eigen::Matrix2d positive_metric =
      Eigen::Vector2d(2.0, 1.0).asDiagonal();
  const Eigen::Vector2d oblique_gradient(1.0, 1.0);
  const NeoProblem oblique_problem = dense_problem(
      positive_hessian, positive_metric, oblique_gradient);
  CiahWorkspace oblique(oblique_problem, Eigen::VectorXd::Unit(2, 0));
  oblique.next(oblique_gradient);
  const auto oblique_ray = oblique.origin_ray(oblique_gradient);
  require(oblique_ray && oblique_ray->residual_norm > 0.5,
      "origin-ray fixture did not retain a full-space residual");
  check_origin_ray_residual(*oblique_ray, positive_hessian,
      positive_metric, oblique_gradient);
}

void check_origin_ray_macro_energy() {
  const Eigen::MatrixXd hessian = Eigen::MatrixXd::Constant(1, 1, -20.0 / 3.0);
  const Eigen::MatrixXd metric = Eigen::MatrixXd::Identity(1, 1);
  const Eigen::VectorXd gradient = Eigen::VectorXd::Ones(1);
  const NeoProblem problem = dense_problem(hessian, metric, gradient);
  CiahWorkspace workspace(problem, gradient);
  workspace.next(gradient);
  const auto ray = workspace.origin_ray(gradient);
  require(ray.has_value(), "origin AH ray missing in macro energy fixture");
  const double max_component = ray->step.cwiseAbs().maxCoeff();
  const double alpha_scale = ray->reference_component > 0.0
      ? 1.0 / ray->reference_component
      : std::numeric_limits<double>::infinity();
  const double factor = std::min(alpha_scale, 0.03 / max_component);
  const double clipped = factor * ray->step[0];
  require(clipped < 0.0 && std::abs(clipped) <= 0.03 + 1.0e-15,
      "origin ray clipping did not give a bounded descending trial");
  const auto energy = [](double x) {
    return x - (10.0 / 3.0) * x * x - 900.0 * x * x * x;
  };
  for (int j = 0; j < 6; ++j) {
    const double fraction = std::ldexp(1.0, -j);
    require(energy(0.03 * fraction) > 0.0 &&
        energy(clipped * fraction) < 0.0,
        "origin ray did not repair the uphill cumulative macro ray");
  }
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
  int metric_calls = 0, norm_calls = 0;
  double last_eigenvalue = 0.0;
  CiahOptions options;
  options.residual_preconditioner =
      [&](const Eigen::VectorXd& residual, double eigenvalue) {
        ++preconditioner_calls;
        last_eigenvalue = eigenvalue;
        return Eigen::VectorXd(residual);
      };
  const NeoProblem problem = dense_problem(
      hessian, metric, old_gradient, &actions, &metric_calls, &norm_calls);
  CiahWorkspace workspace(problem, Eigen::VectorXd::Unit(2, 0), options);
  const CiahResult first = workspace.next(old_gradient);
  require(first.residual_norm < 1.0e-12 && actions == 1,
          "one-dimensional AH residual should vanish with one action");
  const CiahResult settled = workspace.next(old_gradient);
  require(settled.converged && settled.new_hvp_actions == 0,
          "unchanged AH root did not satisfy eigenvalue-change convergence");
  workspace.reset_rhs();
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
  require(metric_calls == 2 && norm_calls == 0,
          "expansion diagnostics evaluated an extra metric or inverse norm");
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

/** @brief An early projected AH step may increase the full model gradient. */
void check_fixed_rhs_forcing_and_reprojection() {
  constexpr double epsilon = 0.01;
  constexpr double coupling = 4.0;
  Eigen::Matrix2d hessian;
  hessian << 1.0, coupling,
             coupling, coupling * coupling + 1.0;
  const Eigen::Matrix2d metric = Eigen::Matrix2d::Identity();
  const Eigen::Vector2d gradient(epsilon, 0.0);
  int actions = 0, metric_calls = 0, norm_calls = 0;
  const NeoProblem problem = dense_problem(hessian, metric, gradient,
      &actions, &metric_calls, &norm_calls);
  CiahWorkspace workspace(problem, Eigen::VectorXd::Unit(2, 0));
  const CiahResult early = workspace.next(gradient);
  check_images_and_residual(hessian, metric, gradient, early);
  const Eigen::VectorXd changed_gradient = gradient + early.hessian_step;
  require(hessian.determinant() > 0.0 && early.eigenvalue < 0.0 &&
          early.step.norm() < 0.03 &&
          gradient.dot(early.step) +
              0.5 * early.step.dot(early.hessian_step) < 0.0 &&
          changed_gradient.norm() > 3.0 * gradient.norm(),
          "SPD projected AH fixture did not increase the model residual");
  require(!early.resolved(gradient.norm(), 1.0 / 3.0) &&
          early.resolved(gradient.norm(), 10.0) &&
          !early.resolved(-1.0, 1.0 / 3.0) &&
          !early.resolved(gradient.norm(), -1.0),
          "AH forcing did not use the full unscaled residual");

  // Consuming this deliberately unresolved step changes the AH RHS. The old
  // pending residual must not cause a Hessian action before a new Ritz solve.
  workspace.reset_rhs();
  const CiahResult reprojected = workspace.next(changed_gradient);
  require(reprojected.new_hvp_actions == 0 && actions == 1 &&
          metric_calls == 1 && norm_calls == 0 &&
          reprojected.basis_size == early.basis_size && !reprojected.converged,
          "RHS refresh expanded an obsolete residual or reused convergence history");
  check_images_and_residual(hessian, metric, changed_gradient, reprojected);
  const Eigen::VectorXd obsolete_rhs_residual = reprojected.reference_component *
      (gradient + hessian * reprojected.step -
       reprojected.eigenvalue * metric * reprojected.step);
  require((reprojected.residual - obsolete_rhs_residual).norm() > epsilon,
          "RHS refresh retained the old AH right-hand side");
  const CiahResult enriched = workspace.next(changed_gradient);
  require(enriched.new_hvp_actions == 1 && actions == 2 &&
          metric_calls == 2 && norm_calls == 0 &&
          enriched.resolved(changed_gradient.norm(), 1.0 / 3.0),
          "fixed refreshed RHS did not enrich to a resolved full-space AH root");
  check_images_and_residual(hessian, metric, changed_gradient, enriched);
  check_full_space_root(hessian, metric, changed_gradient, enriched);
  const CiahResult settled = workspace.next(changed_gradient);
  require(settled.converged && settled.new_hvp_actions == 0 && actions == 2,
          "fixed refreshed RHS repeated an action after full-space convergence");

  // Waiting on the original fixed RHS gives the independent dense AH oracle.
  int original_actions = 0;
  const NeoProblem original = dense_problem(
      hessian, metric, gradient, &original_actions);
  CiahWorkspace fixed(original, Eigen::VectorXd::Unit(2, 0));
  const CiahResult unresolved = fixed.next(gradient);
  require(!unresolved.resolved(gradient.norm(), 1.0 / 3.0),
          "fixed-RHS oracle fixture was already resolved in one direction");
  const CiahResult resolved = fixed.next(gradient);
  require(resolved.resolved(gradient.norm(), 1.0 / 3.0) &&
          original_actions == 2,
          "fixed-RHS AH enrichment failed the relative forcing condition");
  check_full_space_root(hessian, metric, gradient, resolved);

  CiahResult small_rhs;
  small_rhs.residual_norm = 0.0;
  require(!small_rhs.resolved(0.0, 1.0 / 3.0),
          "zero reference weight passed unscaled AH forcing");
  small_rhs.reference_component = 0.2;
  require(small_rhs.resolved(0.0, 1.0 / 3.0),
          "exact AH equation failed zero-RHS forcing");
  small_rhs.residual_norm = 1.0e-12;
  require(!small_rhs.resolved(0.0, 1.0 / 3.0) &&
          small_rhs.resolved(1.0e-9, 1.0 / 3.0),
          "small-RHS forcing lost its reference-component scaling");
}

/** @brief An RHS-specific dependent expansion must not exhaust future RHSs. */
void check_rhs_reset_clears_expansion_exhaustion() {
  const Eigen::Matrix2d hessian = Eigen::Vector2d(2.0, 3.0).asDiagonal();
  const Eigen::Matrix2d metric = Eigen::Matrix2d::Identity();
  const Eigen::Vector2d old_gradient(0.2, 0.5);
  const Eigen::Vector2d new_gradient(0.3, 0.4);
  int actions = 0;
  const NeoProblem problem = dense_problem(hessian, metric, old_gradient, &actions);
  bool dependent = true;
  CiahOptions options;
  options.residual_preconditioner =
      [&](const Eigen::VectorXd&, double) {
        Eigen::VectorXd direction = Eigen::VectorXd::Zero(2);
        direction[dependent ? 0 : 1] = 1.0;
        return direction;
      };
  options.residual_direction = [](const Eigen::VectorXd&) {
    return Eigen::VectorXd::Unit(2, 0);
  };
  CiahWorkspace workspace(problem, Eigen::VectorXd::Unit(2, 0), options);
  const CiahResult first = workspace.next(old_gradient);
  require(first.residual_norm > 0.1 && actions == 1,
      "dependent-expansion fixture has no unresolved complement");
  const CiahResult exhausted = workspace.next(old_gradient);
  require(exhausted.exhaustion && exhausted.new_hvp_actions == 0 && actions == 1,
      "dependent expansion did not mark its original RHS exhausted");
  dependent = false;
  workspace.reset_rhs();
  const CiahResult reprojected = workspace.next(new_gradient);
  require(!reprojected.exhaustion && reprojected.new_hvp_actions == 0 &&
          actions == 1,
      "new RHS inherited an old dependent-expansion exhaustion");
  check_images_and_residual(hessian, metric, new_gradient, reprojected);
  const CiahResult enriched = workspace.next(new_gradient);
  require(enriched.new_hvp_actions == 1 && actions == 2 &&
          enriched.basis_size == 2,
      "new RHS could not enrich after old dependent expansion");
  check_full_space_root(hessian, metric, new_gradient, enriched);
  workspace.reset_rhs();
  const CiahResult capped = workspace.next(old_gradient);
  require(capped.exhaustion && capped.new_hvp_actions == 0 && actions == 2,
      "RHS reset bypassed the hard basis-size cap");
}

/** @brief Incremental projections agree with full rebuilding after every admission. */
void check_incremental_projection() {
  constexpr int n = 7;
  Eigen::MatrixXd h = Eigen::MatrixXd::Zero(n,n);
  Eigen::MatrixXd g = Eigen::MatrixXd::Zero(n,n);
  Eigen::VectorXd gradient(n);
  for (int i=0;i<n;++i) {
    gradient[i] = 0.3/(i+1);
    for (int j=0;j<n;++j) {
      h(i,j) = 0.2/(i+j+1) + (i==j ? i+1.0 : 0.0);
      g(i,j) = 0.1/(i+j+1) + (i==j ? 0.2*(i+1) : 0.0);
    }
  }
  int actions = 0;
  std::vector<Eigen::VectorXd> directions;
  CiahOptions options;
  options.max_cycle = n;
  options.observe_images = [&](const Eigen::VectorXd& v,const Eigen::VectorXd&) {
    directions.push_back(v);
  };
  const auto problem = dense_problem(h,g,gradient,&actions);
  CiahWorkspace workspace(problem,gradient,options);
  for (int iteration=0;iteration<n+3;++iteration) {
    // Refresh the RHS even after full-space exhaustion: frozen projections
    // must remain reusable without either a Hessian or metric evaluation.
    Eigen::VectorXd rhs = gradient;
    rhs[iteration%n] += 0.01*(iteration+1);
    if (iteration) workspace.reset_rhs();
    const auto reprojected = workspace.next(rhs);
    if (iteration)
      require(reprojected.new_hvp_actions == 0,
          "changed RHS caused an action before cached-space reprojection");
    const auto result = workspace.next(rhs);
    const int m = directions.size();
    Eigen::MatrixXd q(n,m);
    for (int i=0;i<m;++i) q.col(i)=directions[i];
    Eigen::MatrixXd a=Eigen::MatrixXd::Zero(m+1,m+1);
    a.bottomRightCorner(m,m)=q.transpose()*h*q;
    a.block(1,0,m,1)=q.transpose()*rhs;
    a.block(0,1,1,m)=rhs.transpose()*q;
    Eigen::MatrixXd metric=Eigen::MatrixXd::Identity(m+1,m+1);
    metric.bottomRightCorner(m,m)=q.transpose()*g*q;
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> reference(a,metric);
    require(reference.info()==Eigen::Success,"rebuilt CIAH reference failed");
    int selected=0;
    while (selected<=m && std::abs(reference.eigenvectors()(0,selected))<=0.1)
      ++selected;
    require(selected<=m,"rebuilt CIAH reference has no selected root");
    const Eigen::VectorXd step=q*reference.eigenvectors().col(selected).tail(m)/
        reference.eigenvectors()(0,selected);
    require_close(result.step,step,3.0e-11,"incremental AH projection changed a Ritz step");
    check_images_and_residual(h,g,rhs,result);
    require(actions==m,"incremental AH projection repeated a Hessian action");
  }
  // An actually asymmetric action must still fail the two-sided certificate.
  Eigen::Matrix2d nonsymmetric;
  nonsymmetric << 2.0,0.7,0.0,3.0;
  const auto bad=dense_problem(nonsymmetric,Eigen::Matrix2d::Identity(),
      Eigen::Vector2d(1.0,1.0));
  CiahWorkspace invalid(bad,Eigen::VectorXd::Unit(2,0));
  invalid.next(Eigen::Vector2d(1.0,1.0));
  bool rejected=false;
  try { invalid.next(Eigen::Vector2d(1.0,1.0)); }
  catch (const std::runtime_error&) { rejected=true; }
  require(rejected,"incremental AH projection concealed an asymmetric action");
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
    check_restart_phase();
    workspace.next(gradient);
  } catch (const std::runtime_error&) {
    rejected = true;
  }
  require(rejected, "CIAH accepted a rank-deficient projected metric");
}

/** @brief Keyframe corrections are not part of the cached frozen H*p image. */
void check_accumulated_microstep_images() {
  Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();
  hessian.diagonal() << 2.0, 4.0, 7.0;
  const Eigen::Matrix3d metric = Eigen::Matrix3d::Identity();
  const Eigen::Vector3d initial(0.3, -0.2, 0.1);
  const Eigen::Vector3d correction(0.01, 0.02, -0.03);
  int actions = 0;
  const NeoProblem problem = dense_problem(hessian, metric, initial, &actions);
  CiahWorkspace workspace(problem, initial);
  Eigen::VectorXd displacement = Eigen::VectorXd::Zero(3);
  Eigen::VectorXd image = displacement;
  Eigen::VectorXd predicted = initial;
  Eigen::VectorXd prefix, prefix_image;
  for (int micro = 0; micro < 3; ++micro) {
    if (micro) workspace.reset_rhs();
    const auto candidate = workspace.next(predicted);
    const double scale = micro == 0 ? 0.5 : 0.75;
    displacement += scale * candidate.step;
    image += scale * candidate.hessian_step;
    predicted += scale * candidate.hessian_step;
    if (micro == 0) predicted += correction;
    if (micro == 1) {
      prefix = displacement;
      prefix_image = image;
    }
  }
  const int previous_actions = actions;
  const Eigen::VectorXd frozen_residual = initial + image;
  require_close(image, hessian * displacement, 2.0e-12,
                "accumulated microstep image does not equal H*p");
  require_close(frozen_residual, initial + hessian * displacement, 2.0e-12,
                "macro frozen residual is inconsistent");
  require_close(predicted - frozen_residual, correction, 2.0e-12,
                "keyframe correction was confused with the frozen residual");
  require_close(prefix_image, hessian * prefix, 2.0e-12,
                "retained prefix image is inconsistent");
  require(actions == previous_actions && actions == workspace.basis_size(),
          "macro residual diagnostics repeated a cached Hessian action");
}

/** @brief Energy anchoring must follow the same keyframe as the refreshed RHS. */
void check_keyframe_segment_minimum() {
  using xmvb::vb::CiahEnergyModel;
  for (const double cubic : {0.0, 1e-15, 0.4}) {
    const double minimum = 0.3;
    const auto value = [&](double t) {
      const double x = t-minimum;
      return -10.0+x*x+cubic*x*x*x;
    };
    const auto derivative = [&](double t) {
      const double x = t-minimum;
      return 2*x+3*cubic*x*x;
    };
    // The Hermite data, not the deliberately wrong frozen curvature,
    // must determine the finite-step correction.
    CiahEnergyModel model(value(0), Eigen::VectorXd::Constant(1, derivative(0)));
    model.add(Eigen::VectorXd::Ones(1), Eigen::VectorXd::Constant(1, -7.0));
    const auto trial = model.minimum(value(1), Eigen::VectorXd::Constant(1, derivative(1)));
    require(std::abs(trial.fraction-minimum) < 1e-12,
            "keyframe Hermite minimum depends on stale curvature");
    require(std::abs(trial.energy-value(minimum)) < 1e-12,
            "keyframe Hermite energy is incorrect");
  }
  for (const double slope : {-1.0, 1.0}) {
    CiahEnergyModel model(-10, Eigen::VectorXd::Constant(1, slope));
    model.add(Eigen::VectorXd::Ones(1), Eigen::VectorXd::Zero(1));
    const auto trial = model.minimum(-10+slope, Eigen::VectorXd::Constant(1, slope));
    require(trial.fraction == (slope < 0 ? 1.0 : 0.0),
            "linear segment chose the wrong endpoint");
  }
}

void check_keyframe_energy_model() {
  using xmvb::vb::CiahEnergyModel;
  const Eigen::Vector2d g0(1.0, -0.5);
  Eigen::Matrix2d h;
  h << 2.0, 0.2, 0.2, 3.0;
  CiahEnergyModel model(-10.0, g0);
  const Eigen::Vector2d first(-0.1, 0.05);
  model.add(first, h * first);
  require(std::abs(model.energy() -
      (-10.0 + g0.dot(first) + 0.5 * first.dot(h * first))) < 1e-13,
      "CIAH origin energy model is inconsistent");
  require(model.anchor_energy() == -10.0 &&
      model.descent_mismatch(-9.99) && model.descent_mismatch(-10.0) &&
      !model.descent_mismatch(-10.08),
      "CIAH failed to distinguish true descent from an energy-sign mismatch");
  const Eigen::Vector2d gkf(0.4, -0.1);
  model.refresh(-10.08, gkf);
  require(model.displacement().norm() == 0.0,
      "CIAH keyframe did not reset its interval displacement");
  require(std::abs(model.energy() + 10.08) < 1e-13,
      "CIAH keyframe kept an obsolete energy anchor");
  const Eigen::Vector2d second(-0.03, 0.02);
  model.add(second, h * second);
  require_close(model.displacement(), second, 1e-13,
      "CIAH keyframe interval includes the preceding macro displacement");
  require(std::abs(model.energy() -
      (-10.08 + gkf.dot(second) + 0.5 * second.dot(h * second))) < 1e-13,
      "CIAH model includes displacement preceding the keyframe");
  require(model.anchor_energy() == -10.08 &&
      model.descent_mismatch(-10.075),
      "CIAH compared a keyframe segment with the obsolete macro energy");
  model.remove(second, h * second);
  require(std::abs(model.energy() + 10.08) < 1e-13,
      "CIAH rejected microstep was not removed from the energy model");
  require(!model.descent_mismatch(-10.07),
      "CIAH invented a descent prediction for an empty keyframe segment");
}

/** @brief A converged AH root must not create another old-residual action. */
void check_converged_root_reuse() {
  Eigen::Matrix3d h;
  h << 2.0, 0.1, 0.0, 0.1, 4.0, 0.2, 0.0, 0.2, 6.0;
  const Eigen::Matrix3d m = Eigen::Matrix3d::Identity();
  const Eigen::Vector3d g(1.0, 0.0, 0.0);
  int actions = 0;
  const auto problem = dense_problem(h, m, g, &actions);
  CiahOptions options;
  options.tolerance = 1.0;  // Deliberately converge with a nonzero residual.
  CiahWorkspace workspace(problem, g, options);
  workspace.next(g);
  const auto converged = workspace.next(g);
  require(converged.converged, "AH reuse test did not reach a converged root");
  require(converged.residual_norm > options.lindep,
      "AH reuse test needs a nonzero residual before full-space exhaustion");
  const int before = actions;
  const auto repeated = workspace.next(g);
  require(repeated.new_hvp_actions == 0 && actions == before,
      "converged AH root expanded an obsolete residual");
}

/** @brief Backtracking retains actual descent and matching immutable images. */
void check_keyframe_checkpoint() {
  xmvb::vb::CiahKeyframes keyframes(-10.0);
  require(!keyframes.available(),"empty keyframe history is available");
  const Eigen::Vector2d p(.1,-.2), hp(.3,-.8), last(.02,-.01);
  CiahResult ritz;
  ritz.residual_norm=.02;
  ritz.eigenvalue=-.004;
  keyframes.save(-9.9,p,hp,last,ritz);
  require(!keyframes.available(),"uphill keyframe became a descent checkpoint");
  keyframes.save(-10.1,p,hp,last,ritz);
  keyframes.save(-10.05,2*p,2*hp,2*last,ritz);
  const auto point=keyframes.take();
  require(point.energy == -10.1 && point.residual_norm == ritz.residual_norm &&
      point.eigenvalue == ritz.eigenvalue,"keyframe checkpoint metadata is inconsistent");
  require_close(point.displacement,p,1e-14,"keyframe displacement was overwritten");
  require_close(point.image,hp,1e-14,"keyframe HVP image was overwritten");
  require_close(point.last_step,last,1e-14,"keyframe restart direction was overwritten");
  require(!keyframes.available(),"consumed keyframe can be retried indefinitely");
}

/** @brief Verify symmetry, exact coarse inverse and shifted signed models. */
void check_coarse_balancing_inverse() {
  for (const double first : {2.0,-2.0}) {
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(5,5);
    h.diagonal() << first,3,5,4,6;
    h(0,3)=h(3,0)=0.7;
    h(1,4)=h(4,1)=0.5;
    h(2,3)=h(3,2)=0.2;
    h(2,0)=h(0,2)=0.1;
    Eigen::MatrixXd g = Eigen::VectorXd::LinSpaced(5,1,2).asDiagonal();
    g(0,1)=g(1,0)=0.1;
    xmvb::vb::CoarsePreconditioner coarse(5,4,1.0e-14);
    Eigen::MatrixXd v(5,4);
    v.col(0)=Eigen::VectorXd::Unit(5,0);
    v.col(1)=Eigen::VectorXd::Unit(5,1);
    v.col(2)=Eigen::VectorXd::Unit(5,3);
    v.col(3)=Eigen::VectorXd::Unit(5,4);
    coarse.observe(v.col(0),h*v.col(0),g*v.col(0));
    coarse.observe(v.col(0),h*v.col(0),g*v.col(0));
    require(coarse.size() == 1,"coarse dependent sample was admitted");
    const Eigen::VectorXd mixed=v.col(0)+v.col(1);
    coarse.observe(mixed,h*mixed,g*mixed);
    for (int i=2;i<4;++i) coarse.observe(v.col(i),h*v.col(i),g*v.col(i));
    const auto sampled = coarse.orbital_spectrum(5,
        [&](const Eigen::VectorXd& rhs) -> Eigen::VectorXd { return h.fullPivLu().solve(rhs); });
    require(sampled.dimension == 4 && sampled.positive_curvature == (first > 0),
        "sampled orbital curvature rank or sign is incorrect");
    if (first > 0)
      require(std::abs(sampled.preconditioned_min-1.0)<2.0e-12 &&
          std::abs(sampled.preconditioned_max-1.0)<2.0e-12,
          "exact orbital inverse is not the identity in the sampled spectrum");
    const Eigen::MatrixXd base=Eigen::VectorXd::LinSpaced(5,0.2,0.8).asDiagonal();
    for (const double omega : {0.0,-0.2}) {
      const Eigen::MatrixXd l=h-omega*g;
      const Eigen::MatrixXd e=v.transpose()*l*v;
      const Eigen::MatrixXd coarse_inverse=v*e.inverse()*v.transpose();
      const Eigen::MatrixXd left=Eigen::MatrixXd::Identity(5,5)-coarse_inverse*l;
      const Eigen::MatrixXd expected=coarse_inverse+left*base*left.transpose();
      Eigen::MatrixXd actual(5,5);
      int calls=0;
      for (int i=0;i<5;++i) actual.col(i)=coarse.apply(Eigen::VectorXd::Unit(5,i),omega,
          [&](const Eigen::VectorXd& r) -> Eigen::VectorXd { ++calls;return base*r; });
      require(calls == 5,"coarse inverse repeated the base preconditioner");
      require((actual-expected).norm()<2.0e-12,"coarse inverse differs from dense reference");
      require((actual-actual.transpose()).norm()<2.0e-12,"coarse inverse lost symmetry");
      require((actual*l*v-v).norm()<2.0e-12,"coarse sampled inverse is not exact");
    }
  }
}

/** @brief Cached-image rank tests use Gram eigenvalues, not singular values. */
void check_coarse_sample_rank() {
  xmvb::vb::CoarsePreconditioner coarse(3,3,1.0e-14);
  const Eigen::Matrix3d h=Eigen::Vector3d(2.0,3.0,4.0).asDiagonal();
  const Eigen::Vector3d first=Eigen::Vector3d::Unit(0);
  coarse.observe(first,h*first,first);
  const Eigen::Vector3d unresolved=first+1.0e-8*Eigen::Vector3d::Unit(1);
  coarse.observe(unresolved,h*unresolved,unresolved);
  require(coarse.size()==1,
      "coarse inverse admitted a sample with an unresolved Gram eigenvalue");
  const Eigen::Vector3d resolved=first+1.0e-3*Eigen::Vector3d::Unit(1);
  coarse.observe(resolved,h*resolved,resolved);
  require(coarse.size()==2,"coarse inverse rejected a resolved independent sample");
  xmvb::vb::CoarsePreconditioner chained(3,3,1.0e-14);
  chained.observe(first,h*first,first);
  const Eigen::Vector3d second=first+1.0e-4*Eigen::Vector3d::Unit(1);
  chained.observe(second,h*second,second);
  require(chained.size()==2,"coarse inverse rejected a resolved sample pair");
  const Eigen::Vector3d third=Eigen::Vector3d::Unit(1)+
      1.0e-4*Eigen::Vector3d::Unit(2);
  chained.observe(third,h*third,third);
  require(chained.size()==2,
      "coarse inverse certified local rank but missed global sample dependence");
}

void check_projected_dual_norm() {
  const Eigen::Matrix3d h = Eigen::Vector3d(2.0,3.0,4.0).asDiagonal();
  Eigen::Matrix3d g;
  g << 2.0,0.4,0.1, 0.4,3.0,0.2, 0.1,0.2,4.0;
  const Eigen::Vector3d gradient(0.1,0.2,0.3);
  int actions = 0;
  const auto problem = dense_problem(h,g,gradient,&actions);
  const Eigen::Vector3d seed(1.0,2.0,0.0);
  CiahWorkspace workspace(problem,Eigen::VectorXd(seed));
  require(workspace.projected_dual_norm(gradient) == 0.0,
      "empty CIAH space captures a nonzero covector norm");
  workspace.next(gradient);
  const Eigen::Vector3d covector(0.4,-0.1,0.8);
  const double expected = std::abs(seed.dot(covector))/std::sqrt(seed.dot(g*seed));
  require(std::abs(workspace.projected_dual_norm(covector)-expected)<2.0e-14,
      "CIAH projected norm differs from independent metric projection");
  require(actions == 1,"CIAH projected norm evaluated an extra Hessian image");
}

void check_sampled_schur_spectrum() {
  Eigen::Matrix4d h=Eigen::Vector4d(4.0,3.0,2.0,5.0).asDiagonal();
  h(0,2)=h(2,0)=1.0;
  h(1,3)=h(3,1)=0.5;
  const Eigen::Matrix4d g=Eigen::Vector4d(2.0,3.0,1.5,0.5).asDiagonal();
  for (bool positive : {true,false}) {
    if (!positive) h(3,3)=-5.0;
    xmvb::vb::CoarsePreconditioner coarse(4,4,1.0e-14);
    for (int j=0;j<4;++j) {
      const Eigen::Vector4d v=Eigen::Vector4d::Unit(j);
      coarse.observe(v,h*v,g*v);
    }
    const auto result=coarse.schur_spectrum(2);
    require(result.orbital_dimension==2 && result.structure_dimension==2,
        "sampled Schur split has incorrect dimensions");
    require(result.structure_positive==positive,
        "sampled Schur diagnostic misclassified structure curvature");
    if (!positive) continue;
    require(std::abs(result.orbital_min-1.0)<2.0e-14 &&
        std::abs(result.orbital_max-2.0)<2.0e-14,
        "sampled orbital spectrum ignores its physical metric");
    require(std::abs(result.structure_min-4.0/3.0)<2.0e-14 &&
        std::abs(result.structure_max-10.0)<2.0e-14,
        "sampled structure spectrum ignores its physical metric");
    require(std::abs(result.relaxed_min-2.95/3.0)<2.0e-14 &&
        std::abs(result.relaxed_max-1.75)<2.0e-14 &&
        std::abs(result.maximum_fraction_removed-0.125)<2.0e-14,
        "sampled Schur spectrum differs from independent elimination");
  }
}

/** @brief Cached-image diagnostics match an independent full-space spectrum. */
void check_inverse_spectrum() {
  Eigen::Matrix3d h;
  h << 2.0,0.2,0.1,0.2,5.0,0.3,0.1,0.3,9.0;
  const Eigen::Matrix3d g=Eigen::Vector3d(0.2,1.0,3.0).asDiagonal();
  const Eigen::Vector3d gradient(0.2,0.3,0.7);
  int actions=0;
  const auto problem=dense_problem(h,g,gradient,&actions);
  CiahWorkspace workspace(problem,gradient);
  for (int i=0;i<3;++i) workspace.next(gradient);
  require(workspace.basis_size()==3,"inverse spectrum fixture did not span the space");
  const int before=actions;
  const double shift=-0.07;
  const Eigen::Matrix3d l=h-shift*g;
  const Eigen::LDLT<Eigen::Matrix3d> factor(l);
  const auto exact=workspace.inverse_spectrum(shift,[&](const Eigen::VectorXd& r) {
    return Eigen::VectorXd(factor.solve(r));
  });
  require(exact.positive_curvature && std::abs(exact.inverse_min-1.0)<2.0e-11 &&
      std::abs(exact.inverse_max-1.0)<2.0e-11 && exact.inverse_skew<2.0e-12,
      "exact inverse diagnostic is not the identity spectrum");
  const auto identity=workspace.inverse_spectrum(shift,
      [](const Eigen::VectorXd& r) { return r; });
  const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> reference(l);
  const Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix3d> physical(l,g);
  require(std::abs(identity.inverse_min-reference.eigenvalues()[0])<2.0e-10 &&
      std::abs(identity.inverse_max-reference.eigenvalues()[2])<2.0e-10 &&
      std::abs(identity.curvature_min-physical.eigenvalues()[0])<2.0e-10 &&
      std::abs(identity.curvature_max-physical.eigenvalues()[2])<2.0e-10,
      "cached inverse spectrum differs from independent dense reference");
  const auto negative=workspace.inverse_spectrum(shift,
      [&](const Eigen::VectorXd& r) { return Eigen::VectorXd(-factor.solve(r)); });
  require(std::abs(negative.inverse_min+1.0)<2.0e-11 &&
      std::abs(negative.inverse_max+1.0)<2.0e-11,
      "diagnostic concealed an indefinite inverse");
  int calls=0;
  const auto indefinite=workspace.inverse_spectrum(100.0,[&](const Eigen::VectorXd& r) {
    ++calls; return r;
  });
  require(!indefinite.positive_curvature && !calls,
      "inverse diagnostic treated indefinite curvature as SPD");
  require(actions==before,"inverse spectrum diagnostic performed a Hessian action");
}

/** @brief Norm control includes the complement without solving the physical metric. */
void check_coordinate_norms() {
  Eigen::Matrix3d hessian;
  hessian << 1.0, 0.8, 0.3, 0.8, 2.0, 0.4, 0.3, 0.4, 3.0;
  Eigen::Matrix3d metric;
  metric << 2.0, 0.2, 0.1, 0.2, 1.0, 0.0, 0.1, 0.0, 0.5;
  const Eigen::Vector3d gradient(1.0, 0.0, 0.0);
  int norms = 0;
  const auto problem = dense_problem(hessian, metric, gradient, nullptr, nullptr, &norms);
  CiahWorkspace workspace(problem);
  auto result = workspace.next(gradient);
  require(norms == 0 && !result.converged,
      "CIAH norm control forced a metric inverse or falsely converged");
  require(workspace.projected_dual_norm(result.residual) < 1e-12 &&
      result.residual.norm() > 0.1,
      "coordinate-norm test does not exercise a residual outside the AH subspace");
  check_images_and_residual(hessian, metric, gradient, result);
  auto next = workspace.next(gradient);
  require(next.basis_size == 2, "full residual did not expand the AH space");
  check_images_and_residual(hessian, metric, gradient, next);
  const Eigen::Vector3d refreshed = gradient + Eigen::Vector3d(.02, -.03, .01);
  workspace.reset_rhs();
  next = workspace.next(refreshed);
  require(next.new_hvp_actions == 0 && next.basis_size == 2,
      "coordinate-norm RHS refresh did not reproject cached images first");
  check_images_and_residual(hessian, metric, refreshed, next);
  next = workspace.next(refreshed);
  check_full_space_root(hessian, metric, refreshed, next);
  require(norms == 0, "Ritz or keyframe RHS update evaluated the physical dual norm");
}

}  // namespace

int main() {
  try {
    check_anisotropic_spd();
    check_indefinite_root_selection();
    check_cached_origin_ray();
    check_origin_ray_macro_energy();
    check_gradient_refresh_and_preconditioner();
    check_fixed_rhs_forcing_and_reprojection();
    check_rhs_reset_clears_expansion_exhaustion();
    check_incremental_projection();
    check_metric_validation();
    check_accumulated_microstep_images();
    check_keyframe_energy_model();
    check_keyframe_segment_minimum();
    check_converged_root_reuse();
    check_coarse_balancing_inverse();
    check_keyframe_checkpoint();
    check_coarse_sample_rank();
    check_projected_dual_norm();
    check_sampled_schur_spectrum();
    check_inverse_spectrum();
    check_restart_phase();
    check_coordinate_norms();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
