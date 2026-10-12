#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

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
                            MetricSolveStats* stats = nullptr) {
  if (stats) *stats = {};
  const auto apply = [&](const Eigen::VectorXd& v) {
    if (stats) ++stats->actions;
    return action(v);
  };
  const auto certify = [&](const Eigen::VectorXd& v) {
    if (stats) ++stats->certifications;
    Eigen::VectorXd residual = rhs - apply(v);
    if (stats) stats->relative_residual = residual.norm() / rhs.norm();
    return residual;
  };
  if (!rhs.allFinite()) throw std::invalid_argument("non-finite metric RHS");
  Eigen::VectorXd x = Eigen::VectorXd::Zero(rhs.size());
  if (rhs.size() == 0 || rhs.isZero()) return x;
  const double rhs_norm = rhs.norm();
  const double target = std::sqrt(std::numeric_limits<double>::epsilon()) * rhs_norm;
  Eigen::VectorXd r = rhs, z = inverse(r), p = z;
  double rz = r.dot(z), cycle_norm = rhs_norm;
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

}  // namespace xmvb::core
