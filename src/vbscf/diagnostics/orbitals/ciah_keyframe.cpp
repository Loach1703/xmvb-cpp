#include "vbscf/diagnostics/orbitals/ciah_keyframe.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "vbscf/diagnostics/performance/neo_trace.hpp"

namespace xmvb::vb::diagnostics {

bool ciah_keyframe_audit_enabled() {
  const char* flag = std::getenv("XMVB_NEO_KEYFRAME_AUDIT");
  return neo_event_trace_enabled() && flag &&
      (std::strcmp(flag, "transport") == 0 || std::strcmp(flag, "1") == 0);
}

bool ciah_model_audit_enabled() {
  const char* flag = std::getenv("XMVB_NEO_KEYFRAME_AUDIT");
  return neo_event_trace_enabled() && flag && std::strcmp(flag, "model") == 0;
}

void audit_ciah_model(
    const Eigen::VectorXd& anchor, const CiahKeyframeProbe& origin,
    const Eigen::VectorXd& step, const Eigen::VectorXd& frozen_image,
    const CiahKeyframeProbe& endpoint, int orbital_size,
    const std::function<double(const Eigen::VectorXd&)>& norm,
    const std::function<CiahKeyframeProbe(const Eigen::VectorXd&)>& evaluate) {
  const double slope = origin.gradient.dot(step);
  Eigen::VectorXd local_image;
  for (const double epsilon : {1.0e-3, 3.0e-4}) {
    const auto plus = evaluate(anchor + epsilon * step);
    const auto minus = evaluate(anchor - epsilon * step);
    local_image = (plus.gradient - minus.gradient) / (2.0 * epsilon);
    neo_trace("ciah_model_curvature_audit",
        "anchor_norm=%.17g,epsilon=%.17g,derivative_error=%.17g,local_image=%.17g,frozen_image=%.17g,image_change=%.17g,local_curvature=%.17g,frozen_curvature=%.17g",
        anchor.norm(), epsilon, std::abs((plus.energy-minus.energy)/(2*epsilon)-slope),
        norm(local_image), norm(frozen_image), norm(local_image-frozen_image),
        step.dot(local_image), step.dot(frozen_image));
  }
  for (const double fraction : {1.0, 0.5, 0.25, 0.125}) {
    const auto value = fraction == 1.0 ? endpoint : evaluate(anchor + fraction * step);
    const double actual = value.energy - origin.energy;
    const double frozen = fraction*slope + 0.5*fraction*fraction*step.dot(frozen_image);
    const double local = fraction*slope + 0.5*fraction*fraction*step.dot(local_image);
    const Eigen::VectorXd error = value.gradient-origin.gradient-fraction*frozen_image;
    neo_trace("ciah_model_scale_audit",
        "anchor_norm=%.17g,fraction=%.17g,actual_change=%.17g,frozen_change=%.17g,local_change=%.17g,energy_error=%.17g,local_energy_error=%.17g,gradient_error=%.17g,local_gradient_error=%.17g,orbital_error_euclidean=%.17g,structure_error_euclidean=%.17g",
        anchor.norm(), fraction, actual, frozen, local, actual-frozen, actual-local,
        norm(error), norm(value.gradient-origin.gradient-fraction*local_image),
        error.head(orbital_size).norm(), error.tail(error.size()-orbital_size).norm());
  }
  Eigen::VectorXd orbital = step, structure = step;
  orbital.tail(step.size()-orbital_size).setZero();
  structure.head(orbital_size).setZero();
  const auto orbital_point = evaluate(anchor+orbital);
  const auto structure_point = evaluate(anchor+structure);
  neo_trace("ciah_model_block_audit",
      "anchor_norm=%.17g,orbital_change=%.17g,structure_change=%.17g,coupled_change=%.17g,mixed_change=%.17g",
      anchor.norm(), orbital_point.energy-origin.energy, structure_point.energy-origin.energy,
      endpoint.energy-origin.energy,
      endpoint.energy-orbital_point.energy-structure_point.energy+origin.energy);
}

void audit_ciah_keyframe(
    const Eigen::VectorXd& displacement,
    const Eigen::VectorXd& hessian_displacement,
    const Eigen::VectorXd& actual_gradient,
    int orbital_size,
    const std::function<CiahKeyframeProbe(const Eigen::VectorXd&)>& evaluate) {
  for (const double epsilon : {1.0e-3, 3.0e-4}) {
    const auto plus = evaluate(epsilon * displacement);
    const auto minus = evaluate(-epsilon * displacement);
    const Eigen::VectorXd error = (plus.gradient - minus.gradient) / (2.0 * epsilon) -
        hessian_displacement;
    neo_trace("ciah_transported_hvp_audit",
        "epsilon=%.17g,full_error=%.17g,orbital_error=%.17g,structure_error=%.17g,image_norm=%.17g",
        epsilon, error.norm(), error.head(orbital_size).norm(),
        error.tail(error.size() - orbital_size).norm(), hessian_displacement.norm());
    const auto finite_plus = evaluate((1.0 + epsilon) * displacement);
    const auto finite_minus = evaluate((1.0 - epsilon) * displacement);
    const double derivative = actual_gradient.dot(displacement);
    const double finite = (finite_plus.energy - finite_minus.energy) / (2.0 * epsilon);
    neo_trace("ciah_transported_gradient_audit", "epsilon=%.17g,analytic=%.17g,finite=%.17g,error=%.17g",
        epsilon, derivative, finite, std::abs(derivative - finite));
  }
}

}  // namespace xmvb::vb::diagnostics
