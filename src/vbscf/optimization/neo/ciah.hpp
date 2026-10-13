#pragma once

#include <functional>
#include <optional>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "vbscf/optimization/neo/problem.hpp"

namespace xmvb::vb {

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

  /** @brief Incremental displacement since the last accepted keyframe. */
  const Eigen::VectorXd& displacement() const noexcept { return step_; }

private:
  double energy_;
  Eigen::VectorXd gradient_, step_, image_;
};

/** Controls one bounded, matrix-free CIAH Davidson workspace. */
struct CiahOptions {
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
  /** @brief Euclidean norm of the complete residual, including outside the trial space. */
  double residual_norm = 0.0;
  /** @brief Euclidean norm of gradient + H * step, before clipping. */
  double model_gradient_norm = 0.0;
  bool converged = false;
  int basis_size = 0;
  /** Hessian actions made by this call; zero or one. */
  int new_hvp_actions = 0;
  /** No further subspace expansion is available. */
  bool exhaustion = false;
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
 * The NeoProblem must outlive this workspace. `next` may receive a refreshed
 * gradient after a nonlinear microstep; all admitted directions and their
 * Hessian/metric images remain attached to the same immutable problem.
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
      const NeoProblem& problem,
      std::optional<Eigen::VectorXd> initial_direction = std::nullopt,
      CiahOptions options = {});

  explicit CiahWorkspace(const NeoProblem& problem, CiahOptions options)
      : CiahWorkspace(problem, std::nullopt, std::move(options)) {}

  /** Append the previous residual direction, then return the next AH root. */
  CiahResult next(const Eigen::VectorXd& current_gradient);

  int basis_size() const noexcept {
    return static_cast<int>(basis_.size());
  }

private:
  Eigen::VectorXd precondition(
      const Eigen::VectorXd& residual, double eigenvalue) const;
  bool append(Eigen::VectorXd direction,
      const Eigen::VectorXd* residual = nullptr,
      double residual_norm = 0.0, double shift = 0.0);

  const NeoProblem& problem_;
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
  bool exhausted_ = false;
};

}  // namespace xmvb::vb
