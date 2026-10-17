#pragma once

#include <functional>
#include <Eigen/Core>
#include "linalg/actions.hpp"

namespace xmvb::mcscf {

/** @brief Covector norm supplied by the physical model, not by the solver. */
using DualNorm = std::function<double(const Eigen::VectorXd&)>;

} // namespace xmvb::mcscf
