#pragma once

#include <functional>
#include <Eigen/Core>

namespace xmvb::mcscf {

/** @brief Matrix-free vector-to-vector action in a caller-defined coordinate space. */
using VectorAction = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;

/** @brief Covector norm supplied by the physical model, not by the solver. */
using DualNorm = std::function<double(const Eigen::VectorXd&)>;

} // namespace xmvb::mcscf
