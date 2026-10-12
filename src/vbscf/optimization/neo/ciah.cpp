#include "vbscf/optimization/neo/ciah.hpp"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <limits>
#include <stdexcept>
#include <utility>

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#include "vbscf/optimization/neo/symmetry.hpp"
#include "vbscf/diagnostics/performance/neo_trace.hpp"

namespace xmvb::vb {
namespace {

Eigen::MatrixXd columns(const std::vector<Eigen::VectorXd>& vectors) {
  Eigen::MatrixXd result(vectors.front().size(), vectors.size());
  for (Eigen::Index j = 0; j < result.cols(); ++j) {
    result.col(j) = vectors[static_cast<std::size_t>(j)];
  }
  return result;
}

void check_projected_symmetry(
    const Eigen::MatrixXd& basis,
    const Eigen::MatrixXd& images,
    const Eigen::MatrixXd& action,
    const char* name,
    double relative_action_error = 0.0) {
  const Eigen::MatrixXd skew = action - action.transpose();
  const double roundoff = neo_detail::projection_skew_roundoff(
      basis, images, action);
  const double tolerance = 256.0 * std::numeric_limits<double>::epsilon() *
      static_cast<double>(action.rows()) * std::max(1.0, action.norm()) +
      roundoff + relative_action_error * basis.norm() * images.norm();
  if (!std::isfinite(skew.norm()) || skew.norm() > tolerance) {
    throw std::runtime_error(name);
  }
}

}  // namespace

CiahWorkspace::CiahWorkspace(
    const NeoProblem& problem,
    std::optional<Eigen::VectorXd> initial_direction,
    CiahOptions options)
    : problem_(problem),
      options_(std::move(options)),
      initial_direction_(std::move(initial_direction)) {
  if (!(options_.tolerance > 0.0) || !std::isfinite(options_.tolerance) ||
      !(options_.lindep > 0.0) || !std::isfinite(options_.lindep) ||
      options_.lindep >= 1.0 || options_.max_cycle < 1 ||
      options_.operator_relative_accuracy < 0.0 ||
      !std::isfinite(options_.operator_relative_accuracy)) {
    throw std::invalid_argument("invalid CIAH Davidson options");
  }
  if (initial_direction_ &&
      (initial_direction_->size() != problem_.size() ||
       !initial_direction_->allFinite())) {
    throw std::invalid_argument("invalid CIAH initial direction");
  }
}

Eigen::VectorXd CiahWorkspace::precondition(
    const Eigen::VectorXd& residual, double eigenvalue) const {
  const bool trace = diagnostics::neo_event_trace_enabled();
  const auto start = trace ? std::chrono::steady_clock::now() :
      std::chrono::steady_clock::time_point{};
  Eigen::VectorXd direction = options_.residual_preconditioner
      ? options_.residual_preconditioner(residual, eigenvalue)
      : problem_.has_preconditioner()
          ? problem_.apply_preconditioner(residual)
          : residual;
  if (direction.size() != problem_.size() || !direction.allFinite()) {
    throw std::runtime_error("CIAH preconditioner returned an invalid vector");
  }
  if (trace) diagnostics::neo_trace("ciah_precondition_cost",
      "basis=%d,seconds=%.17g", basis_size(),
      std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
  return direction;
}

bool CiahWorkspace::append(Eigen::VectorXd direction) {
  if (direction.size() != problem_.size() || !direction.allFinite()) {
    throw std::invalid_argument("invalid CIAH Davidson direction");
  }
  const double initial_norm = direction.stableNorm();
  if (!(initial_norm > 0.0)) return false;
  const bool trace = diagnostics::neo_event_trace_enabled();
  const auto start = trace ? std::chrono::steady_clock::now() :
      std::chrono::steady_clock::time_point{};
  // Euclidean two-pass orthogonalization avoids another Hessian action.
  for (int pass = 0; pass < 2; ++pass) {
    for (const Eigen::VectorXd& vector : basis_) {
      direction.noalias() -= vector.dot(direction) * vector;
    }
  }
  const double norm = direction.stableNorm();
  if (!(norm > options_.lindep * initial_norm) ||
      !(norm > 64.0 * std::numeric_limits<double>::epsilon() *
          initial_norm)) {
    return false;
  }
  direction /= norm;
  const auto metric_start = trace ? std::chrono::steady_clock::now() :
      std::chrono::steady_clock::time_point{};
  // Evaluate each image once, only after the direction has been admitted.
  Eigen::VectorXd metric_image = problem_.apply_metric(direction);
  const auto hessian_start = trace ? std::chrono::steady_clock::now() :
      std::chrono::steady_clock::time_point{};
  Eigen::VectorXd hessian_image = problem_.apply_hessian(direction);
  if (options_.observe_images) options_.observe_images(direction,metric_image);
  if (trace) diagnostics::neo_trace("ciah_expansion_cost",
      "basis=%d,remaining_fraction=%.17g,orthogonalize_seconds=%.17g,metric_seconds=%.17g,hvp_seconds=%.17g",
      basis_size()+1, norm/initial_norm,
      std::chrono::duration<double>(metric_start-start).count(),
      std::chrono::duration<double>(hessian_start-metric_start).count(),
      std::chrono::duration<double>(std::chrono::steady_clock::now()-hessian_start).count());
  basis_.push_back(std::move(direction));
  metric_images_.push_back(std::move(metric_image));
  hessian_images_.push_back(std::move(hessian_image));
  return true;
}

CiahResult CiahWorkspace::next(const Eigen::VectorXd& current_gradient) {
  if (current_gradient.size() != problem_.size() ||
      !current_gradient.allFinite()) {
    throw std::invalid_argument("invalid CIAH current gradient");
  }
  CiahResult result;
  if (basis_.empty()) {
    Eigen::VectorXd seed = initial_direction_
        ? *initial_direction_
        : precondition(-current_gradient, 0.0);
    if (!append(std::move(seed))) {
      // A stationary gradient still needs a probe to form an AH Ritz root.
      Eigen::VectorXd probe = Eigen::VectorXd::Ones(problem_.size());
      if (!append(std::move(probe))) {
        throw std::runtime_error("CIAH could not construct an initial basis");
      }
    }
    result.new_hvp_actions = 1;
    initial_direction_.reset();
  } else if (pending_residual_ && !exhausted_) {
    if (basis_size() >= std::min<int>(options_.max_cycle, problem_.size())) {
      exhausted_ = true;
    } else {
      Eigen::VectorXd direction = precondition(
          *pending_residual_, pending_eigenvalue_);
      if (!append(std::move(direction)) &&
          !append(*pending_residual_)) {
        exhausted_ = true;
      } else {
        result.new_hvp_actions = 1;
      }
    }
  }
  pending_residual_.reset();

  const bool trace = diagnostics::neo_event_trace_enabled();
  const auto projection_start = trace ? std::chrono::steady_clock::now() :
      std::chrono::steady_clock::time_point{};
  const Eigen::MatrixXd q = columns(basis_);
  const Eigen::MatrixXd hq = columns(hessian_images_);
  const Eigen::MatrixXd gq = columns(metric_images_);
  const Eigen::MatrixXd projected_hessian_action = q.transpose() * hq;
  const Eigen::MatrixXd projected_metric_action = q.transpose() * gq;
  check_projected_symmetry(
      q, hq, projected_hessian_action,
      "CIAH Hessian action is not symmetric",
      options_.operator_relative_accuracy);
  check_projected_symmetry(
      q, gq, projected_metric_action,
      "CIAH metric action is not symmetric");

  const Eigen::MatrixXd projected_hessian = 0.5 *
      (projected_hessian_action + projected_hessian_action.transpose());
  const Eigen::MatrixXd projected_metric = 0.5 *
      (projected_metric_action + projected_metric_action.transpose());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> metric_spectrum(
      projected_metric, Eigen::EigenvaluesOnly);
  if (metric_spectrum.info() != Eigen::Success ||
      !(metric_spectrum.eigenvalues()[0] >
        options_.lindep * metric_spectrum.eigenvalues().tail(1)[0])) {
    throw std::runtime_error("CIAH projected metric is not positive definite or is rank deficient");
  }

  const Eigen::Index m = q.cols();
  Eigen::MatrixXd augmented_hessian = Eigen::MatrixXd::Zero(m + 1, m + 1);
  const Eigen::VectorXd projected_gradient = q.transpose() * current_gradient;
  augmented_hessian.block(0, 1, 1, m) = projected_gradient.transpose();
  augmented_hessian.block(1, 0, m, 1) = projected_gradient;
  augmented_hessian.bottomRightCorner(m, m) = projected_hessian;
  Eigen::MatrixXd augmented_metric = Eigen::MatrixXd::Identity(m + 1, m + 1);
  augmented_metric.bottomRightCorner(m, m) = projected_metric;

  const auto spectrum_start = trace ? std::chrono::steady_clock::now() :
      std::chrono::steady_clock::time_point{};
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(
      augmented_hessian, augmented_metric);
  if (spectrum.info() != Eigen::Success) {
    throw std::runtime_error("CIAH augmented eigensolver failed");
  }
  Eigen::Index selected = -1;
  for (Eigen::Index i = 0; i < spectrum.eigenvalues().size(); ++i) {
    if (std::abs(spectrum.eigenvectors()(0, i)) > 0.1) {
      selected = i;
      break;
    }
  }
  if (selected < 0) {
    throw std::runtime_error("CIAH found no augmented root with sufficient reference weight");
  }

  const double alpha = spectrum.eigenvectors()(0, selected);
  const Eigen::VectorXd coefficients =
      spectrum.eigenvectors().col(selected).tail(m) / alpha;
  result.eigenvalue = spectrum.eigenvalues()[selected];
  result.reference_weight = std::abs(alpha);
  result.step = q * coefficients;
  result.hessian_step = hq * coefficients;
  result.metric_step = gq * coefficients;
  result.residual = alpha *
      (current_gradient + result.hessian_step -
       result.eigenvalue * result.metric_step);
  const auto norm_start = trace ? std::chrono::steady_clock::now() :
      std::chrono::steady_clock::time_point{};
  result.residual_norm = problem_.dual_norm(result.residual);
  // g + Hp = r/alpha + omega Gp. Its complete dual norm needs no second
  // inverse: ||r/alpha||_G^-1^2 + 2 omega r^T p/alpha + omega^2 p^T Gp.
  // Retain the cross term rather than assuming exact Galerkin orthogonality.
  const double residual_length = result.residual_norm / alpha;
  const double squared_model_norm = residual_length * residual_length +
      2.0 * result.eigenvalue * result.residual.dot(result.step) / alpha +
      result.eigenvalue * result.eigenvalue * result.step.dot(result.metric_step);
  if (!(squared_model_norm >= 0.0) || !std::isfinite(squared_model_norm)) {
    throw std::runtime_error("CIAH model gradient norm is invalid");
  }
  result.model_gradient_norm = std::sqrt(squared_model_norm);
  if (trace) diagnostics::neo_trace("ciah_ritz_cost",
      "basis=%d,projection_seconds=%.17g,spectrum_image_seconds=%.17g,norm_seconds=%.17g,gradient_euclidean=%.17g,residual_euclidean=%.17g,galerkin_residual=%.17g",
      basis_size(), std::chrono::duration<double>(spectrum_start - projection_start).count(),
      std::chrono::duration<double>(norm_start - spectrum_start).count(),
      std::chrono::duration<double>(std::chrono::steady_clock::now() - norm_start).count(),
      current_gradient.norm(), result.residual.norm(), (q.transpose()*result.residual).norm());
  result.converged = previous_eigenvalue_ &&
      std::abs(result.eigenvalue - *previous_eigenvalue_) < options_.tolerance &&
      result.residual_norm < std::sqrt(options_.tolerance);
  previous_eigenvalue_ = result.eigenvalue;
  result.basis_size = basis_size();
  result.exhaustion = exhausted_ ||
      basis_size() >= std::min<int>(options_.max_cycle, problem_.size());
  // davidson_cc does not append a residual after yielding a converged root.
  // Reproject a changed RHS in the retained space first; that new root can
  // supply a new residual if it is no longer converged.
  if (!result.converged && !result.exhaustion &&
      result.residual_norm >= options_.lindep) {
    pending_residual_ = result.residual;
    pending_eigenvalue_ = result.eigenvalue;
  }
  return result;
}

}  // namespace xmvb::vb
