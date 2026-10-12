#pragma once

#include <functional>
#include <Eigen/Core>

namespace xmvb::vb {

/** @brief Symmetric balancing inverse on an exact sampled coupled space.
 * Images and metric images must belong to the same frozen accepted point.
 * No operator evaluation or cross-point history is retained by this class.
 */
class CoarsePreconditioner {
public:
  CoarsePreconditioner(int dimension, int capacity, double lindep, double action_accuracy = 0.0);
  /** @brief Orthogonalize a sample and both images with identical coefficients. */
  void observe(const Eigen::VectorXd& direction, const Eigen::VectorXd& image,
      const Eigen::VectorXd& metric_image);
  /** @brief Apply the signed balancing inverse with one base inverse call. */
  Eigen::VectorXd apply(const Eigen::VectorXd& residual, double eigenvalue,
      const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& inverse) const;
  int size() const noexcept { return size_; }

private:
  int size_ = 0;
  double lindep_, action_accuracy_;
  Eigen::MatrixXd basis_, images_, metrics_;
};

}  // namespace xmvb::vb
