#include "mcscf/optimization/ciah/solver.hpp"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#include "mcscf/optimization/trace.hpp"
#include "core/timing.hpp"

namespace xmvb::mcscf {
namespace {

void check_projected_symmetry(
    Eigen::Index dimension,
    double basis_norm, double images_norm,
    const Eigen::MatrixXd& magnitude,
    const Eigen::MatrixXd& action,
    const char* name,
    double relative_action_error = 0.0) {
  const Eigen::MatrixXd skew = action - action.transpose();
  // The same dot-product certificate as projection_skew_roundoff, with
  // positive dot sums cached when each row/column was first constructed.
  const double u = 0.5 * std::numeric_limits<double>::epsilon();
  const double nu = static_cast<double>(dimension) * u;
  const double gamma = nu / (1.0 - nu);
  const Eigen::MatrixXd error =
      (gamma / (1.0 - gamma)) * (magnitude + magnitude.transpose()) +
      u * (action.cwiseAbs() + action.transpose().cwiseAbs());
  const double roundoff = error.stableNorm();
  const double tolerance = 256.0 * std::numeric_limits<double>::epsilon() *
      static_cast<double>(action.rows()) * std::max(1.0, action.norm()) +
      roundoff + relative_action_error * basis_norm * images_norm;
  if (!std::isfinite(skew.norm()) || skew.norm() > tolerance) {
    std::ostringstream message;
    message.precision(17);
    message << name << ": skew=" << skew.norm()
            << ", tolerance=" << tolerance << ", dot_roundoff=" << roundoff
            << ", action_norm=" << action.norm()
            << ", basis_norm=" << basis_norm << ", images_norm=" << images_norm
            << ", dimension=" << dimension << ", subspace=" << action.rows();
    throw std::runtime_error(message.str());
  }
}

}  // namespace

bool CiahResult::resolved(double rhs_norm, double relative_accuracy) const noexcept {
  const double alpha = std::abs(reference_component);
  return std::isfinite(residual_norm) && residual_norm >= 0.0 &&
      std::isfinite(alpha) && alpha > 0.0 &&
      std::isfinite(rhs_norm) && rhs_norm >= 0.0 &&
      std::isfinite(relative_accuracy) && relative_accuracy >= 0.0 &&
      residual_norm <= alpha * relative_accuracy * rhs_norm;
}

CiahWorkspace::CiahWorkspace(
    const QuadraticProblem& problem,
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
  const Eigen::Index capacity = std::min<Eigen::Index>(options_.max_cycle, problem_.size());
  hessian_projection_ = Eigen::MatrixXd::Zero(capacity, capacity);
  metric_projection_ = hessian_projection_;
  hessian_magnitude_ = hessian_projection_;
  metric_magnitude_ = hessian_projection_;
}

void CiahWorkspace::reset_rhs() noexcept {
  pending_residual_.reset();
  pending_eigenvalue_ = 0.0;
  pending_residual_norm_ = 0.0;
  previous_eigenvalue_.reset();
  // An expansion can be dependent for one RHS but useful for the next.
  // A hard basis-size cap is checked again by next().
  exhausted_ = false;
}

Eigen::VectorXd CiahWorkspace::precondition(
    const Eigen::VectorXd& residual, double eigenvalue) const {
  xmvb::profiling::ScopedPhase phase(xmvb::profiling::Phase::NeoPreconditioner);
  const bool trace = diagnostics::ciah_event_trace_enabled();
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
  if (trace) diagnostics::ciah_trace("ciah_precondition_cost",
      "basis=%d,seconds=%.17g", basis_size(),
      std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
  return direction;
}

bool CiahWorkspace::append(Eigen::VectorXd direction,
    const Eigen::VectorXd* residual, double residual_norm, double shift) {
  xmvb::profiling::ScopedPhase basis_phase(xmvb::profiling::Phase::NeoBasis);
  if (direction.size() != problem_.size() || !direction.allFinite()) {
    throw std::invalid_argument("invalid CIAH Davidson direction");
  }
  const double initial_norm = direction.stableNorm();
  if (!(initial_norm > 0.0)) return false;
  const bool trace = diagnostics::ciah_event_trace_enabled();
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
  basis_phase.stop();
  xmvb::profiling::ScopedPhase images_phase(xmvb::profiling::Phase::CiahImages, true);
  const auto metric_start = trace ? std::chrono::steady_clock::now() :
      std::chrono::steady_clock::time_point{};
  // Evaluate each image once, only after the direction has been admitted.
  Eigen::VectorXd metric_image, hessian_image;
  if (!options_.apply_images) metric_image = problem_.apply_metric(direction);
  const auto hessian_start = trace ? std::chrono::steady_clock::now() :
      std::chrono::steady_clock::time_point{};
  if (options_.apply_images) {
    auto images = options_.apply_images(direction);
    hessian_image = std::move(images.first);
    metric_image = std::move(images.second);
    if (hessian_image.size() != direction.size() ||
        metric_image.size() != direction.size() ||
        !hessian_image.allFinite() || !metric_image.allFinite())
      throw std::runtime_error("invalid joint CIAH operator images");
  } else {
    hessian_image = problem_.apply_hessian(direction);
  }
  images_phase.stop();
  if (diagnostics::ciah_trace_enabled() && residual) {
    // Audit BEFORE observe_images learns this direction. Afterwards coarse
    // exactness would make the inverse test tautological on its new samples.
    // This is diagnostic-only: reuse Lq, add one inverse/metric action, never
    // use the error to change a candidate, tolerance, or stopping decision.
    const auto audit_start = std::chrono::steady_clock::now();
    const Eigen::VectorXd shifted_image = hessian_image - shift * metric_image;
    const Eigen::VectorXd error = precondition(shifted_image, shift) - direction;
    Eigen::VectorXd complement = error;
    for (int pass = 0; pass < 2; ++pass)
      for (const auto& vector : basis_)
        complement.noalias() -= vector.dot(complement) * vector;
    const Eigen::VectorXd error_metric = problem_.apply_metric(error);
    const double squared = error.dot(error_metric);
    const double q_squared = direction.dot(metric_image);
    diagnostics::ciah_trace("ciah_unlearned_inverse",
        "basis=%d,shift=%.17g,relative_physical_error=%.17g,error_euclidean=%.17g,complement_euclidean=%.17g,curvature=%.17g,seconds=%.17g",
        basis_size()+1,shift,std::sqrt(std::max(0.0,squared)/q_squared),
        error.norm(),complement.norm(),direction.dot(shifted_image),
        std::chrono::duration<double>(std::chrono::steady_clock::now()-audit_start).count());
  }
  if (options_.observe_images) options_.observe_images(direction,metric_image);
  if (trace && residual && residual_norm > 0.0) {
    // Reuse the progress norm and the admitted direction's exact images.
    // This measures enrichment, not convergence or a global condition number.
    const double metric_squared = direction.dot(metric_image);
    const double work = residual->dot(direction);
    const double curvature = direction.dot(hessian_image) - shift * metric_squared;
    diagnostics::ciah_trace("ciah_expansion_geometry",
        "basis=%d,residual_coordinate_norm=%.17g,euclidean_cosine=%.17g,metric_squared=%.17g,shifted_curvature=%.17g,physical_rayleigh=%.17g,remaining_fraction_estimate=%.17g",
        basis_size()+1,residual_norm,
        work/(residual->norm()*direction.norm()),metric_squared,
        curvature,curvature/metric_squared,norm/initial_norm);
  }
  if (trace) diagnostics::ciah_trace("ciah_expansion_cost",
      "basis=%d,remaining_fraction=%.17g,orthogonalize_seconds=%.17g,metric_seconds=%.17g,hvp_seconds=%.17g",
      basis_size()+1, norm/initial_norm,
      std::chrono::duration<double>(metric_start-start).count(),
      std::chrono::duration<double>(hessian_start-metric_start).count(),
      std::chrono::duration<double>(std::chrono::steady_clock::now()-hessian_start).count());
  basis_.push_back(std::move(direction));
  metric_images_.push_back(std::move(metric_image));
  hessian_images_.push_back(std::move(hessian_image));
  const auto projection_start = trace ? std::chrono::steady_clock::now() :
      std::chrono::steady_clock::time_point{};
  const Eigen::Index j = basis_size() - 1;
  xmvb::profiling::ScopedPhase projection_phase(xmvb::profiling::Phase::NeoProjectionFactor);
  for (Eigen::Index i = 0; i <= j; ++i) {
    // Keep both orientations: symmetrizing an assumed-symmetric triangle
    // would conceal inconsistent matrix-free images.
    hessian_projection_(i,j) = basis_[i].dot(hessian_images_[j]);
    hessian_projection_(j,i) = basis_[j].dot(hessian_images_[i]);
    metric_projection_(i,j) = basis_[i].dot(metric_images_[j]);
    metric_projection_(j,i) = basis_[j].dot(metric_images_[i]);
    hessian_magnitude_(i,j) = basis_[i].cwiseAbs().dot(hessian_images_[j].cwiseAbs());
    hessian_magnitude_(j,i) = basis_[j].cwiseAbs().dot(hessian_images_[i].cwiseAbs());
    metric_magnitude_(i,j) = basis_[i].cwiseAbs().dot(metric_images_[j].cwiseAbs());
    metric_magnitude_(j,i) = basis_[j].cwiseAbs().dot(metric_images_[i].cwiseAbs());
  }
  basis_squared_norm_ += basis_[j].squaredNorm();
  hessian_squared_norm_ += hessian_images_[j].squaredNorm();
  metric_squared_norm_ += metric_images_[j].squaredNorm();
  if (trace) diagnostics::ciah_trace("ciah_projection_update_cost",
      "basis=%d,seconds=%.17g", basis_size(),
      std::chrono::duration<double>(std::chrono::steady_clock::now()-projection_start).count());
  return true;
}

double CiahWorkspace::projected_dual_norm(const Eigen::VectorXd& covector) const {
  if (covector.size() != problem_.size() || !covector.allFinite())
    throw std::invalid_argument("invalid CIAH diagnostic covector");
  if (basis_.empty()) return 0.0;
  if (projected_metric_.rows() != basis_size())
    throw std::logic_error("CIAH diagnostic requires a completed Ritz solve");
  Eigen::VectorXd rhs(basis_size());
  for (int i = 0; i < basis_size(); ++i) rhs[i] = basis_[i].dot(covector);
  const Eigen::LLT<Eigen::MatrixXd> factor(projected_metric_);
  if (factor.info() != Eigen::Success)
    throw std::runtime_error("CIAH diagnostic metric factorization failed");
  const Eigen::VectorXd solution = factor.solve(rhs);
  const double squared = rhs.dot(solution);
  if (!(squared >= 0.0) || !std::isfinite(squared))
    throw std::runtime_error("CIAH diagnostic projection norm is invalid");
  return std::sqrt(squared);
}

CiahInverseSpectrum CiahWorkspace::inverse_spectrum(double shift,
    const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& inverse) const {
  if (!std::isfinite(shift) || !inverse)
    throw std::invalid_argument("invalid CIAH inverse spectrum input");
  CiahInverseSpectrum result;
  const int m = basis_size();
  result.dimension = m;
  if (!m) return result;
  const Eigen::MatrixXd raw = hessian_projection_.topLeftCorner(m,m) -
      shift * metric_projection_.topLeftCorner(m,m);
  const Eigen::MatrixXd e = 0.5 * (raw+raw.transpose());
  const Eigen::MatrixXd g = 0.5 *
      (metric_projection_.topLeftCorner(m,m) +
       metric_projection_.topLeftCorner(m,m).transpose());
  const Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> curvature(
      e,g,Eigen::EigenvaluesOnly);
  if (curvature.info() != Eigen::Success)
    throw std::runtime_error("CIAH trial curvature spectrum failed");
  result.curvature_min = curvature.eigenvalues()[0];
  result.curvature_max = curvature.eigenvalues().tail(1)[0];
  result.positive_curvature = result.curvature_min >
      options_.lindep * result.curvature_max;
  if (!result.positive_curvature) return result;
  Eigen::MatrixXd t(m,m);
  for (int j=0;j<m;++j) {
    const Eigen::VectorXd rhs = hessian_images_[j] - shift*metric_images_[j];
    const Eigen::VectorXd image = inverse(rhs);
    if (image.size() != problem_.size() || !image.allFinite())
      throw std::runtime_error("invalid CIAH diagnostic inverse image");
    for (int i=0;i<m;++i)
      t(i,j) = hessian_images_[i].dot(image) - shift*metric_images_[i].dot(image);
  }
  result.inverse_skew = (t-t.transpose()).norm()/std::max(1.0,t.norm());
  t = (0.5*(t+t.transpose())).eval();
  const Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(
      t,e,Eigen::EigenvaluesOnly);
  if (spectrum.info() != Eigen::Success)
    throw std::runtime_error("CIAH trial inverse spectrum failed");
  result.inverse_min = spectrum.eigenvalues()[0];
  result.inverse_max = spectrum.eigenvalues().tail(1)[0];
  return result;
}

std::optional<CiahOriginRay> CiahWorkspace::origin_ray(
    const Eigen::VectorXd& origin_gradient) const {
  if (origin_gradient.size() != problem_.size() ||
      !origin_gradient.allFinite()) {
    throw std::invalid_argument("invalid CIAH origin gradient");
  }
  const Eigen::Index m = basis_size();
  if (!m) return std::nullopt;
  if (projected_metric_.rows() != m || projected_metric_.cols() != m) {
    throw std::logic_error("CIAH origin ray requires a completed Ritz solve");
  }

  // next() has certified the two cached projections and the metric rank.
  const Eigen::MatrixXd raw_hessian = hessian_projection_.topLeftCorner(m, m);
  const Eigen::MatrixXd projected_hessian =
      0.5 * (raw_hessian + raw_hessian.transpose());
  Eigen::VectorXd projected_gradient(m);
  for (Eigen::Index i = 0; i < m; ++i)
    projected_gradient[i] = basis_[i].dot(origin_gradient);
  Eigen::MatrixXd augmented_hessian = Eigen::MatrixXd::Zero(m + 1, m + 1);
  augmented_hessian.block(0, 1, 1, m) = projected_gradient.transpose();
  augmented_hessian.block(1, 0, m, 1) = projected_gradient;
  augmented_hessian.bottomRightCorner(m, m) = projected_hessian;
  Eigen::MatrixXd augmented_metric = Eigen::MatrixXd::Identity(m + 1, m + 1);
  augmented_metric.bottomRightCorner(m, m) = projected_metric_;
  const Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(
      augmented_hessian, augmented_metric);
  if (spectrum.info() != Eigen::Success ||
      !spectrum.eigenvalues().allFinite() ||
      !spectrum.eigenvectors().allFinite()) {
    throw std::runtime_error("CIAH origin augmented eigensolver failed");
  }

  const Eigen::VectorXd coefficients = spectrum.eigenvectors().col(0).tail(m);
  CiahOriginRay ray;
  ray.eigenvalue = spectrum.eigenvalues()[0];
  ray.reference_component = std::abs(spectrum.eigenvectors()(0, 0));
  ray.step = Eigen::VectorXd::Zero(problem_.size());
  ray.hessian_step = ray.step;
  Eigen::VectorXd metric_step = ray.step;
  for (Eigen::Index i = 0; i < m; ++i) {
    ray.step.noalias() += coefficients[i] * basis_[i];
    ray.hessian_step.noalias() += coefficients[i] * hessian_images_[i];
    metric_step.noalias() += coefficients[i] * metric_images_[i];
  }
  const Eigen::VectorXd residual =
      spectrum.eigenvectors()(0, 0) * origin_gradient + ray.hessian_step -
      ray.eigenvalue * metric_step;
  if (!ray.step.allFinite() || !ray.hessian_step.allFinite() ||
      !residual.allFinite())
    throw std::runtime_error("CIAH origin ray is invalid");
  ray.residual_norm = options_.coordinate_norm
      ? options_.coordinate_norm(residual)
      : residual.stableNorm();
  if (!(ray.residual_norm >= 0.0) || !std::isfinite(ray.residual_norm))
    throw std::runtime_error("CIAH origin ray is invalid");
  if (!(ray.step.stableNorm() > 0.0)) return std::nullopt;
  if (origin_gradient.dot(ray.step) > 0.0) {
    ray.step = -ray.step;
    ray.hessian_step = -ray.hessian_step;
  }
  // A zero-slope tail is useful only if it exposes negative curvature.
  if (!(origin_gradient.dot(ray.step) < 0.0 ||
        ray.step.dot(ray.hessian_step) < 0.0)) return std::nullopt;
  return ray;
}

CiahResult CiahWorkspace::next(const Eigen::VectorXd& current_gradient) {
  xmvb::profiling::ScopedPhase next_phase(xmvb::profiling::Phase::CiahNext);
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
      if (!append(std::move(direction), &*pending_residual_,
                  pending_residual_norm_, pending_eigenvalue_) &&
          !append(options_.residual_direction
                      ? options_.residual_direction(*pending_residual_) : *pending_residual_,
                  &*pending_residual_,
                  pending_residual_norm_, pending_eigenvalue_)) {
        exhausted_ = true;
      } else {
        result.new_hvp_actions = 1;
      }
    }
  }
  pending_residual_.reset();

  const bool trace = diagnostics::ciah_event_trace_enabled();
  const auto projection_start = trace ? std::chrono::steady_clock::now() :
      std::chrono::steady_clock::time_point{};
  const Eigen::Index m = basis_size();
  xmvb::profiling::ScopedPhase projection_phase(xmvb::profiling::Phase::NeoProjectedSolve);
  const Eigen::MatrixXd projected_hessian_action = hessian_projection_.topLeftCorner(m,m);
  const Eigen::MatrixXd projected_metric_action = metric_projection_.topLeftCorner(m,m);
  check_projected_symmetry(
      problem_.size(), std::sqrt(basis_squared_norm_), std::sqrt(hessian_squared_norm_),
      hessian_magnitude_.topLeftCorner(m,m), projected_hessian_action,
      "CIAH Hessian action is not symmetric",
      options_.operator_relative_accuracy);
  check_projected_symmetry(
      problem_.size(), std::sqrt(basis_squared_norm_), std::sqrt(metric_squared_norm_),
      metric_magnitude_.topLeftCorner(m,m), projected_metric_action,
      "CIAH metric action is not symmetric");

  const Eigen::MatrixXd projected_hessian = 0.5 *
      (projected_hessian_action + projected_hessian_action.transpose());
  const Eigen::MatrixXd projected_metric = 0.5 *
      (projected_metric_action + projected_metric_action.transpose());
  projected_metric_ = projected_metric;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> metric_spectrum(
      projected_metric, Eigen::EigenvaluesOnly);
  if (metric_spectrum.info() != Eigen::Success ||
      !(metric_spectrum.eigenvalues()[0] >
        options_.lindep * metric_spectrum.eigenvalues().tail(1)[0])) {
    throw std::runtime_error("CIAH projected metric is not positive definite or is rank deficient");
  }

  Eigen::MatrixXd augmented_hessian = Eigen::MatrixXd::Zero(m + 1, m + 1);
  Eigen::VectorXd projected_gradient(m);
  for (Eigen::Index i = 0; i < m; ++i)
    projected_gradient[i] = basis_[i].dot(current_gradient);
  augmented_hessian.block(0, 1, 1, m) = projected_gradient.transpose();
  augmented_hessian.block(1, 0, m, 1) = projected_gradient;
  augmented_hessian.bottomRightCorner(m, m) = projected_hessian;
  Eigen::MatrixXd augmented_metric = Eigen::MatrixXd::Identity(m + 1, m + 1);
  augmented_metric.bottomRightCorner(m, m) = projected_metric;

  const auto spectrum_start = trace ? std::chrono::steady_clock::now() :
      std::chrono::steady_clock::time_point{};
  projection_phase.stop();
  xmvb::profiling::ScopedPhase spectrum_phase(xmvb::profiling::Phase::NeoSpectrum);
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
  result.reference_component = alpha;
  result.step = Eigen::VectorXd::Zero(problem_.size());
  result.hessian_step = result.step;
  result.metric_step = result.step;
  for (Eigen::Index i = 0; i < m; ++i) {
    result.step.noalias() += coefficients[i] * basis_[i];
    result.hessian_step.noalias() += coefficients[i] * hessian_images_[i];
    result.metric_step.noalias() += coefficients[i] * metric_images_[i];
  }
  result.residual = alpha *
      (current_gradient + result.hessian_step -
       result.eigenvalue * result.metric_step);
  spectrum_phase.stop();
  xmvb::profiling::ScopedPhase residual_phase(xmvb::profiling::Phase::NeoProjectionResidual);
  const auto norm_start = trace ? std::chrono::steady_clock::now() :
      std::chrono::steady_clock::time_point{};
  // CIAH uses a coordinate norm to control its inexact microiterations, as
  // PySCF does. The generalized AH metric is still present in the equation.
  // Evaluate the FULL residual, not just its projection into the trial space.
  const auto norm = [&](const Eigen::VectorXd& covector) {
    return options_.coordinate_norm ? options_.coordinate_norm(covector) : covector.norm();
  };
  result.residual_norm = norm(result.residual);
  result.model_gradient_norm = norm(current_gradient + result.hessian_step);
  const bool eigenvalue_settled = previous_eigenvalue_ &&
      std::abs(result.eigenvalue - *previous_eigenvalue_) < options_.tolerance;
  if (trace) {
    Eigen::VectorXd projected_residual(m);
    for (Eigen::Index i = 0; i < m; ++i)
      projected_residual[i] = basis_[i].dot(result.residual);
    diagnostics::ciah_trace("ciah_ritz_cost",
      "basis=%d,projection_seconds=%.17g,spectrum_image_seconds=%.17g,norm_seconds=%.17g,gradient_euclidean=%.17g,residual_euclidean=%.17g,galerkin_residual=%.17g",
      basis_size(), std::chrono::duration<double>(spectrum_start - projection_start).count(),
      std::chrono::duration<double>(norm_start - spectrum_start).count(),
      std::chrono::duration<double>(std::chrono::steady_clock::now() - norm_start).count(),
      current_gradient.norm(), result.residual.norm(), projected_residual.norm());
  }
  result.converged = eigenvalue_settled &&
      result.residual_norm < std::sqrt(options_.tolerance);
  previous_eigenvalue_ = result.eigenvalue;
  result.basis_size = basis_size();
  result.exhaustion = exhausted_ ||
      basis_size() >= std::min<int>(options_.max_cycle, problem_.size());
  // davidson_cc does not append a residual after yielding a converged root.
  // Reproject a changed RHS in the retained space first; that new root can
  // supply a new residual if it is no longer converged.
  // Test backward resolution in the equation's own units, not against an
  // absolute residual amplitude. Rank is then tested on normalized directions.
  const double residual_scale=current_gradient.stableNorm()+
      result.hessian_step.stableNorm()+
      std::abs(result.eigenvalue)*result.metric_step.stableNorm();
  const bool residual_nonzero = result.residual_norm > options_.lindep*residual_scale;
  if (!residual_nonzero && !result.converged) result.exhaustion=true;
  if (!result.converged && !result.exhaustion && residual_nonzero) {
    pending_residual_ = result.residual;
    pending_eigenvalue_ = result.eigenvalue;
    pending_residual_norm_ = result.residual_norm;
  }
  return result;
}


}  // namespace xmvb::mcscf
