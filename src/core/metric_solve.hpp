#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <Eigen/Core>

namespace xmvb::core {

/** @brief Optional accounting; never changes the certified stopping rule. */
struct MetricSolveStats {
  int iterations = 0;
  int actions = 0;
  int certifications = 0;
  double relative_residual = 0.0;
};

/** @brief Matrix-free SPD Riesz solve with a certified true residual.
 * Uses O(n) workspace. A complete CG cycle restarts rather than accepting an
 * inaccurate solution. No approximate inverse is returned on stagnation.
 */
template<class Action, class Inverse>
Eigen::VectorXd solve_metric(const Eigen::VectorXd& rhs,
                            const Action& action, const Inverse& inverse,
                            MetricSolveStats* stats = nullptr,
                            const Eigen::VectorXd* initial = nullptr,
                            Eigen::VectorXd* metric_image = nullptr) {
  if (stats) *stats = {};
  const auto apply = [&](const Eigen::VectorXd& v) {
    if (stats) ++stats->actions;
    return action(v);
  };
  const auto certify = [&](const Eigen::VectorXd& v) {
    if (stats) ++stats->certifications;
    Eigen::VectorXd image = apply(v);
    if (metric_image) *metric_image = image;
    Eigen::VectorXd residual = rhs - image;
    if (stats) stats->relative_residual = residual.norm() / rhs.norm();
    return residual;
  };
  if (!rhs.allFinite()) throw std::invalid_argument("non-finite metric RHS");
  Eigen::VectorXd x = Eigen::VectorXd::Zero(rhs.size());
  if (rhs.size() == 0 || rhs.isZero()) {
    if (metric_image) *metric_image = x;
    return x;
  }
  const double rhs_norm = rhs.norm();
  const double target = std::sqrt(std::numeric_limits<double>::epsilon()) * rhs_norm;
  Eigen::VectorXd r = rhs;
  if (initial) {
    if (initial->size() != rhs.size() || !initial->allFinite())
      throw std::invalid_argument("invalid metric initial solution");
    x = *initial;
    r = certify(x);
    if (r.norm() <= target) return x;
  }
  Eigen::VectorXd z = inverse(r), p = z;
  double rz = r.dot(z), cycle_norm = r.norm();
  Eigen::Index cycle = 0;
  while (true) {
    if (!(rz > 0.0) || !std::isfinite(rz))
      throw std::runtime_error("metric inverse is not positive definite");
    if (stats) ++stats->iterations;
    Eigen::VectorXd ap = apply(p);
    double pap = p.dot(ap);
    if (!(pap > 0.0) || !std::isfinite(pap)) {
      r = certify(x);
      if (r.norm() <= target) return x;
      z = inverse(r);
      p = z;
      ap = apply(p);
      pap = p.dot(ap);
      rz = r.dot(z);
      cycle_norm = r.norm();
      cycle = 0;
      if (!(pap > 0.0) || !std::isfinite(pap) || !(rz > 0.0))
        throw std::runtime_error("metric is singular in tangent space");
    }
    const double alpha = rz / pap;
    if (!std::isfinite(alpha)) throw std::runtime_error("non-finite metric CG step");
    x.noalias() += alpha * p;
    r.noalias() -= alpha * ap;
    if (!r.allFinite()) throw std::runtime_error("non-finite metric residual");
    if (r.norm() <= target || ++cycle == rhs.size()) {
      r = certify(x);
      const double norm = r.norm();
      if (!std::isfinite(norm)) throw std::runtime_error("non-finite true metric residual");
      if (norm <= target) return x;
      if (!(cycle_norm - norm > std::numeric_limits<double>::epsilon() *
            std::max(cycle_norm, rhs_norm)))
        throw std::runtime_error("metric CG stagnated above certified tolerance");
      z = inverse(r);
      p = z;
      rz = r.dot(z);
      cycle_norm = norm;
      cycle = 0;
      continue;
    }
    z = inverse(r);
    const double next_rz = r.dot(z);
    const double beta = next_rz / rz;
    if (!(next_rz > 0.0) || !std::isfinite(beta))
      throw std::runtime_error("invalid metric CG recurrence");
    p = z + beta * p;
    rz = next_rz;
  }
}

/** @brief Fixed-operator Galerkin recycling for certified Riesz solves.
 * Stores at most capacity primal vectors and their certified metric images.
 * The initial solution is Q Q^T rhs with Q^T G Q = I. The balanced SPD inverse
 * is Q Q^T + (I-Q Q^T G) P (I-G Q Q^T), so CG does not relearn the captured
 * metric directions. Every solution passes the unchanged true-residual test.
 * No matrix of size
 * n-by-n is formed, and the workspace must not be shared between operators or
 * concurrent calls.
 */
class MetricSolveWorkspace {
public:
  explicit MetricSolveWorkspace(int capacity) : capacity_(capacity) {
    if (capacity < 0) throw std::invalid_argument("negative metric space capacity");
  }

  int size() const noexcept { return static_cast<int>(basis_.size()); }

  /** @brief Certify a strict dual-norm comparison, including uncaptured RHS.
   * A variational lower bound can reject the comparison after one metric
   * action. It cannot accept it: unresolved comparisons use the unchanged
   * true-residual Riesz solve. No global spectral estimate is assumed.
   */
  template<class Action, class Inverse>
  bool norm_below(const Eigen::VectorXd& rhs, double threshold,
                  const Action& action, const Inverse& inverse) {
    if (threshold <= 0.0) return false;
    if (rhs.isZero()) return true;
    Eigen::VectorXd trial = Eigen::VectorXd::Zero(rhs.size());
    Eigen::VectorXd residual = rhs;
    for (int i = 0; i < size(); ++i) {
      const double weight = basis_[i].dot(rhs);
      trial.noalias() += weight * basis_[i];
      residual.noalias() -= weight * images_[i];
    }
    // The uncaptured residual is retained, not discarded by Galerkin projection.
    trial += inverse(residual);
    const Eigen::VectorXd image = action(trial);
    const double work = rhs.dot(trial), squared = trial.dot(image);
    if (squared > 0.0 && std::abs(work) >= threshold * std::sqrt(squared))
      return false;
    const Eigen::VectorXd solution = solve(rhs, action, inverse);
    return rhs.dot(solution) < threshold * threshold;
  }

  /** @brief Admit an already computed primal/metric pair without a new action. */
  void observe(const Eigen::VectorXd& vector, const Eigen::VectorXd& image) {
    if (size() >= capacity_ || vector.isZero()) return;
    if (image.size() != vector.size() || !vector.allFinite() || !image.allFinite())
      throw std::invalid_argument("invalid recycled metric pair");
    Eigen::VectorXd q = vector, gq = image;
    const double original_squared = q.dot(gq);
    if (!(original_squared > 0.0) || !std::isfinite(original_squared)) return;
    // Reorthogonalize in the physical metric without another operator action.
    for (int pass = 0; pass < 2; ++pass) {
      for (int i = 0; i < size(); ++i) {
        const double weight = basis_[i].dot(gq);
        q.noalias() -= weight * basis_[i];
        gq.noalias() -= weight * images_[i];
      }
    }
    const double squared = q.dot(gq);
    // Numerical dependence only controls reuse, never solve accuracy.
    if (squared > std::numeric_limits<double>::epsilon() * original_squared &&
        std::isfinite(squared)) {
      const double length = std::sqrt(squared);
      basis_.push_back(q / length);
      images_.push_back(gq / length);
    }
  }

  template<class Action, class Inverse>
  Eigen::VectorXd solve(const Eigen::VectorXd& rhs, const Action& action,
                       const Inverse& inverse, MetricSolveStats* stats = nullptr) {
    Eigen::VectorXd initial = Eigen::VectorXd::Zero(rhs.size());
    for (const auto& q : basis_) {
      if (q.size() != rhs.size()) throw std::invalid_argument("metric space size changed");
      initial.noalias() += q.dot(rhs) * q;
    }
    Eigen::VectorXd image;
    const auto balanced_inverse = [&](const Eigen::VectorXd& residual) {
      Eigen::VectorXd projected = residual;
      for (int i = 0; i < size(); ++i)
        projected.noalias() -= basis_[i].dot(residual) * images_[i];
      Eigen::VectorXd z = inverse(projected);
      const Eigen::VectorXd unprojected = z;
      for (int i = 0; i < size(); ++i)
        z.noalias() += (basis_[i].dot(residual) - images_[i].dot(unprojected)) * basis_[i];
      return z;
    };
    Eigen::VectorXd solution = solve_metric(rhs, action, balanced_inverse, stats,
        basis_.empty() ? nullptr : &initial, &image);
    observe(solution, image);
    return solution;
  }

private:
  int capacity_;
  std::vector<Eigen::VectorXd> basis_, images_;
};

}  // namespace xmvb::core
