#pragma once

#include <functional>
#include <optional>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "vbscf/optimization/neo/problem.hpp"

namespace xmvb::vb {

/** Controls one bounded, matrix-free CIAH Davidson workspace. */
struct CiahOptions {
  /** Absolute tolerance in the problem's physical dual norm. */
  double tolerance = 1.0e-12;
  /** Relative threshold for dependent directions and projected metric rank. */
  double lindep = 1.0e-14;
  /** Maximum number of admitted directions (and Hessian actions). */
  int max_cycle = 30;
  /** Declared relative error of the Hessian action for symmetry checks. */
  double operator_relative_accuracy = 0.0;
  /** Optional eigenvalue-shifted residual preconditioner. */
  std::function<Eigen::VectorXd(const Eigen::VectorXd&, double)>
      residual_preconditioner;
};

/** A Ritz candidate, including its exact matrix-free augmented residual. */
struct CiahResult {
  double eigenvalue = 0.0;
  /** Absolute reference component of the generalized-normalized Ritz vector. */
  double reference_weight = 0.0;
  Eigen::VectorXd step;
  Eigen::VectorXd hessian_step;
  Eigen::VectorXd metric_step;
  /** alpha * (gradient + H * step - eigenvalue * G * step). */
  Eigen::VectorXd residual;
  double residual_norm = 0.0;
  bool converged = false;
  int basis_size = 0;
  /** Hessian actions made by this call; zero or one. */
  int new_hvp_actions = 0;
  /** No further subspace expansion is available. */
  bool exhaustion = false;
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
  bool append(Eigen::VectorXd direction);

  const NeoProblem& problem_;
  CiahOptions options_;
  std::optional<Eigen::VectorXd> initial_direction_;
  std::optional<Eigen::VectorXd> pending_residual_;
  std::optional<double> previous_eigenvalue_;
  double pending_eigenvalue_ = 0.0;
  std::vector<Eigen::VectorXd> basis_;
  std::vector<Eigen::VectorXd> hessian_images_;
  std::vector<Eigen::VectorXd> metric_images_;
  bool exhausted_ = false;
};

}  // namespace xmvb::vb
