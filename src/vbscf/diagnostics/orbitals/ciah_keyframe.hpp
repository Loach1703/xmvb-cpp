#pragma once

#include <functional>
#include <Eigen/Core>

namespace xmvb::vb::diagnostics {

struct CiahKeyframeProbe {
  Eigen::VectorXd gradient;
  double energy = 0.0;
};

/** @brief Opt-in derivative audit, never an optimizer acceptance rule. */
bool ciah_keyframe_audit_enabled();

/** @brief Expensive finite-step diagnosis, enabled only by the model audit flag. */
bool ciah_model_audit_enabled();

/** @brief Compare step scaling and local versus frozen curvature at an anchor. */
void audit_ciah_model(
    const Eigen::VectorXd& anchor, const CiahKeyframeProbe& origin,
    const Eigen::VectorXd& step, const Eigen::VectorXd& frozen_image,
    const CiahKeyframeProbe& endpoint, int orbital_size,
    const std::function<double(const Eigen::VectorXd&)>& norm,
    const std::function<CiahKeyframeProbe(const Eigen::VectorXd&)>& evaluate);

/** @brief Opt-in gradient/HVP check at the first keyframe of each macro. */
void audit_ciah_keyframe(
    const Eigen::VectorXd& displacement,
    const Eigen::VectorXd& hessian_displacement,
    const Eigen::VectorXd& actual_gradient,
    int orbital_size,
    const std::function<CiahKeyframeProbe(const Eigen::VectorXd&)>& evaluate);

}  // namespace xmvb::vb::diagnostics
