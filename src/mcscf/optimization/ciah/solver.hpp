#pragma once

#include <functional>
#include <optional>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "mcscf/optimization/problem.hpp"

namespace xmvb::mcscf {

/** @brief Frozen curvature with an energy/gradient anchor refreshed at keyframes.
 *
 * Only displacement since the last keyframe belongs in its quadratic model.
 * The Hessian images remain reusable, as in PySCF's CIAH iteration.
 */
class CiahEnergyModel {
public:
  CiahEnergyModel(double energy, const Eigen::VectorXd& gradient)
      : energy_(energy), gradient_(gradient),
        step_(Eigen::VectorXd::Zero(gradient.size())), image_(step_) {}

  void add(const Eigen::VectorXd& step, const Eigen::VectorXd& image) {
    step_ += step;
    image_ += image;
  }

  void remove(const Eigen::VectorXd& step, const Eigen::VectorXd& image) {
    step_ -= step;
    image_ -= image;
  }

  void refresh(double energy, const Eigen::VectorXd& gradient) {
    energy_ = energy;
    gradient_ = gradient;
    step_.setZero();
    image_.setZero();
  }

  double energy() const {
    return energy_ + gradient_.dot(step_) + 0.5 * step_.dot(image_);
  }

  /** @brief An exact keyframe contradicts the segment's predicted descent.
   * Compare both energies to the last true anchor before refreshing it.
   * This stops extrapolation of an invalid model, not outer convergence.
   */
  bool descent_mismatch(double exact_energy) const {
    return energy() < energy_ && exact_energy >= energy_;
  }

  /** @brief Last true keyframe energy, without subsequent model increments. */
  double anchor_energy() const noexcept { return energy_; }

  /** @brief Incremental displacement since the last accepted keyframe. */
  const Eigen::VectorXd& displacement() const noexcept { return step_; }

private:
  double energy_;
  Eigen::VectorXd gradient_, step_, image_;
};

/** Controls one bounded, matrix-free CIAH Davidson workspace. */
struct CiahOptions {
  /** @brief Joint H/G action when both images share physical intermediates. */
  std::function<std::pair<Eigen::VectorXd, Eigen::VectorXd>(const Eigen::VectorXd&)>
      apply_images;
  /** @brief Local coordinate dual norm; does not solve the physical metric. */
  std::function<double(const Eigen::VectorXd&)> coordinate_norm;
  /** @brief Convert a residual covector to a trial vector if no enrichment remains. */
  std::function<Eigen::VectorXd(const Eigen::VectorXd&)> residual_direction;
  /** Eigenvalue tolerance; the complete coordinate residual uses its square root. */
  double tolerance = 1.0e-12;
  /** Relative threshold for dependent directions and projected metric rank. */
  double lindep = 1.0e-14;
  /** Maximum number of admitted directions (and Hessian actions). */
  int max_cycle = 30;
  /** Declared relative error of the Hessian action for symmetry checks. */
  double operator_relative_accuracy = 0.0;
  /** @brief Observe a direction after its HVP, reusing the computed metric image. */
  std::function<void(const Eigen::VectorXd&,const Eigen::VectorXd&)> observe_images;
  /** Optional eigenvalue-shifted residual preconditioner. */
  std::function<Eigen::VectorXd(const Eigen::VectorXd&, double)>
      residual_preconditioner;
};

/** A Ritz candidate, including its exact matrix-free augmented residual. */
struct CiahResult {
  double eigenvalue = 0.0;
  /** Signed reference component of the generalized-normalized Ritz vector. */
  double reference_component = 0.0;
  Eigen::VectorXd step;
  Eigen::VectorXd hessian_step;
  Eigen::VectorXd metric_step;
  /** alpha * (gradient + H * step - eigenvalue * G * step). */
  Eigen::VectorXd residual;
  /** @brief Coordinate dual norm of the full residual, including outside the trial space. */
  double residual_norm = 0.0;
  /** @brief Coordinate dual norm of gradient + H * step, before clipping. */
  double model_gradient_norm = 0.0;
  bool converged = false;
  int basis_size = 0;
  /** Hessian actions made by this call; zero or one. */
  int new_hvp_actions = 0;
  /** No further subspace expansion is available. */
  bool exhaustion = false;

  /** @brief Inexact solve of the current shifted AH equation.
   * Tests the full residual after removing its reference-component scale.
   * This does not certify a Newton macrostep or decrease of its gradient norm.
   */
  bool resolved(double rhs_norm, double relative_accuracy) const noexcept;
};

/** @brief Direction from the lowest origin-gradient AH root in the retained span.
 *
 * This is an unscaled direction, not a certified Ritz Newton step. The tail
 * is not divided by its reference component, which may be zero for a pure
 * negative-curvature direction. No convergence claim accompanies this ray.
 */
struct CiahOriginRay {
  Eigen::VectorXd step, hessian_step;
  /** Absolute reference weight of the generalized-normalized eigenvector. */
  double reference_component = 0.0;
  double eigenvalue = 0.0;
  /** Coordinate norm of the full augmented residual before clipping. */
  double residual_norm = 0.0;
};

/** @brief Trial-space inverse spectrum, diagnostic only, not a stopping rule. */
struct CiahInverseSpectrum {
  int dimension = 0;
  bool positive_curvature = false;
  double curvature_min = 0.0, curvature_max = 0.0;
  double inverse_min = 0.0, inverse_max = 0.0;
  double inverse_skew = 0.0;
};

/** @brief A true-energy keyframe and its frozen origin-chart images. */
struct CiahCheckpoint {
  double energy = 0.0;
  Eigen::VectorXd displacement, image, last_step;
  double residual_norm = 0.0, eigenvalue = 0.0;
  bool converged = false;
};

/** @brief Retain the lowest evaluated descending keyframe, without new actions. */
class CiahKeyframes {
public:
  explicit CiahKeyframes(double source_energy) : lowest_energy_(source_energy) {}
  void save(double energy, const Eigen::VectorXd& displacement,
      const Eigen::VectorXd& image, const Eigen::VectorXd& last_step,
      const CiahResult& ritz) {
    if (!(energy < lowest_energy_)) return;
    lowest_energy_ = energy;
    checkpoint_ = CiahCheckpoint{energy, displacement, image, last_step,
        ritz.residual_norm, ritz.eigenvalue, ritz.converged};
  }
  bool available() const noexcept { return checkpoint_.has_value(); }
  /** @brief Consume a checkpoint so backtracking cannot retry it indefinitely. */
  CiahCheckpoint take() {
    CiahCheckpoint point = std::move(checkpoint_.value());
    checkpoint_.reset();
    return point;
  }
private:
  double lowest_energy_;
  std::optional<CiahCheckpoint> checkpoint_;
};

/**
 * Davidson iteration for the generalized augmented Hessian of one fixed model.
 *
 * The QuadraticProblem must outlive this workspace. Call `reset_rhs` before supplying
 * a changed gradient to `next`. All admitted directions and their Hessian/metric
 * images remain attached to the same immutable problem.
 */
class CiahWorkspace {
public:
  /** @brief Dual norm captured by the current AH trial space.
   * Diagnostic projection using cached metric images; no operator evaluation.
   */
  double projected_dual_norm(const Eigen::VectorXd& covector) const;
  /** @brief Measure an inverse on cached shifted Hessian images.
   * Computes eigenvalues of (LQ)^T P (LQ) relative to Q^T LQ when the
   * latter is positive definite. Only small matrices and one inverse image
   * are allocated; no Hessian or metric actions are evaluated.
   */
  CiahInverseSpectrum inverse_spectrum(double shift,
      const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& inverse) const;
  explicit CiahWorkspace(
      const QuadraticProblem& problem,
      std::optional<Eigen::VectorXd> initial_direction = std::nullopt,
      CiahOptions options = {});

  explicit CiahWorkspace(const QuadraticProblem& problem, CiahOptions options)
      : CiahWorkspace(problem, std::nullopt, std::move(options)) {}

  /** @brief Mark a changed RHS for a cached-space Ritz solve before expansion.
   * Retains every direction and its H/G images, but discards the old equation's
   * residual and eigenvalue-settling history. With a nonempty retained space,
   * the next call makes no new HVP.
   */
  void reset_rhs() noexcept;

  /** @brief Return the next AH root, enriching on the current RHS residual.
   * Call `reset_rhs` when the supplied gradient changes; that first call only
   * reprojects cached images. Subsequent calls may admit one HVP direction.
   */
  CiahResult next(const Eigen::VectorXd& current_gradient);

  /** @brief Minimize the positive quadratic in the cached span, without HVPs.
   * Returns an unshifted Newton correction only for positive projected
   * curvature and a nonzero descent direction. The residual includes the
   * complement of the span; projected stationarity is not convergence.
   * Call after next(), which validates the frozen H/G projections. The
   * caller owns step feasibility, keyframes and nonlinear acceptance.
   */
  std::optional<CiahResult> projected_newton(const Eigen::VectorXd& gradient);

  /** @brief Re-solve the lowest AH root at an origin gradient using cached data.
   *
   * Requires a completed Ritz solve in the retained span. This does not
   * expand the basis, apply H or G, precondition, or evaluate a physical norm.
   * Returns no ray if that span contains no nonzero descent or negative-
   * curvature direction. Intended for one macro rejection recovery attempt,
   * not as an additional microiteration or convergence criterion.
   */
  std::optional<CiahOriginRay> origin_ray(
      const Eigen::VectorXd& origin_gradient) const;

  int basis_size() const noexcept {
    return static_cast<int>(basis_.size());
  }

private:
  Eigen::VectorXd precondition(
      const Eigen::VectorXd& residual, double eigenvalue) const;
  bool append(Eigen::VectorXd direction,
      const Eigen::VectorXd* residual = nullptr,
      double residual_norm = 0.0, double shift = 0.0);

  const QuadraticProblem& problem_;
  CiahOptions options_;
  std::optional<Eigen::VectorXd> initial_direction_;
  std::optional<Eigen::VectorXd> pending_residual_;
  std::optional<double> previous_eigenvalue_;
  double pending_eigenvalue_ = 0.0;
  /** @brief Cached complete coordinate norm for scalar expansion diagnostics. */
  double pending_residual_norm_ = 0.0;
  std::vector<Eigen::VectorXd> basis_;
  std::vector<Eigen::VectorXd> hessian_images_;
  std::vector<Eigen::VectorXd> metric_images_;
  /** @brief Cached raw projections and positive dot sums for symmetry certification.
   * Only the newly admitted row/column is evaluated; keyframe RHS updates
   * do not modify these frozen-operator projections.
   */
  Eigen::MatrixXd hessian_projection_, metric_projection_;
  Eigen::MatrixXd hessian_magnitude_, metric_magnitude_;
  double basis_squared_norm_ = 0.0;
  double hessian_squared_norm_ = 0.0, metric_squared_norm_ = 0.0;
  Eigen::MatrixXd projected_metric_;
  /** Positive projected factor, reused until a new basis direction arrives. */
  Eigen::MatrixXd newton_lower_;
  int newton_dimension_ = 0;
  bool exhausted_ = false;
};

}  // namespace xmvb::mcscf
