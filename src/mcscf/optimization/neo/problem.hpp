#pragma once

#include <functional>
#include <optional>

#include <Eigen/Core>
#include "mcscf/optimization/actions.hpp"

namespace xmvb::mcscf {

/**
 * @brief Immutable quadratic model used by norm-extended optimization.
 *
 * The Hessian maps a coordinate increment to a gradient covector. The metric
 * defines the physical trust-region norm, and the supplied dual norm measures
 * KKT covectors in the corresponding physical units. An optional SPD
 * preconditioner may accelerate subspace growth, but it is not part of the
 * NEO equations.
 */
class NeoProblem {
public:
  /**
   * @brief Constructs a matrix-free generalized-metric NEO problem.
   *
   * @param gradient Gradient covector at the accepted point.
   * @param apply_hessian Symmetric Hessian action.
   * @param apply_metric Symmetric positive-definite metric action.
   * @param dual_norm Physical dual norm for gradient and KKT covectors.
   * @param apply_preconditioner Optional SPD residual preconditioner.
   * @param hessian_lower_bound Optional certified lower bound on every
   * generalized Hessian eigenvalue. An estimate is not sufficient.
   */
  NeoProblem(
      Eigen::VectorXd gradient,
      VectorAction apply_hessian,
      VectorAction apply_metric,
      DualNorm dual_norm,
      VectorAction apply_preconditioner = {},
      std::optional<double> hessian_lower_bound = std::nullopt);

  /** @brief Number of optimization coordinates. */
  Eigen::Index size() const noexcept { return gradient_.size(); }

  /** @brief Gradient covector at the accepted point. */
  const Eigen::VectorXd& gradient() const noexcept { return gradient_; }

  /** @brief Applies the Hessian and validates its result. */
  Eigen::VectorXd apply_hessian(const Eigen::VectorXd& direction) const;

  /** @brief Applies the coordinate metric and validates its result. */
  Eigen::VectorXd apply_metric(const Eigen::VectorXd& direction) const;

  /** @brief Measures a finite covector in the physical dual metric. */
  double dual_norm(const Eigen::VectorXd& covector) const;

  /** @brief Whether a residual preconditioner was supplied. */
  bool has_preconditioner() const noexcept {
    return static_cast<bool>(apply_preconditioner_);
  }

  /** @brief Applies the optional residual preconditioner. */
  Eigen::VectorXd apply_preconditioner(const Eigen::VectorXd& covector) const;

  /** @brief Certified lower bound on the generalized Hessian spectrum. */
  const std::optional<double>& hessian_lower_bound() const noexcept {
    return hessian_lower_bound_;
  }

private:
  Eigen::VectorXd apply_checked(
      const VectorAction& action,
      const Eigen::VectorXd& vector,
      const char* name) const;

  const Eigen::VectorXd gradient_;
  const VectorAction apply_hessian_;
  const VectorAction apply_metric_;
  const DualNorm dual_norm_;
  const VectorAction apply_preconditioner_;
  const std::optional<double> hessian_lower_bound_;
};

}  // namespace xmvb::mcscf
