#include "vbscf/diagnostics/orbitals/ciah_keyframe.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <Eigen/Cholesky>

#include "vbscf/derivatives/hessian/coupled/structure.hpp"
#include "vbscf/derivatives/hessian/context/accepted_point.hpp"
#include "vbscf/diagnostics/performance/neo_trace.hpp"
#include "vbscf/orbitals/preparation/preparer.hpp"
#include "vbscf/structures/assembly/action.hpp"
#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"
#include "vbscf/structures/orthogonal_ci/metric.hpp"

namespace xmvb::vb::diagnostics {

bool ciah_keyframe_audit_enabled() {
  static const bool enabled = [] {
    const char* flag = std::getenv("XMVB_NEO_KEYFRAME_AUDIT");
    return flag && std::strcmp(flag, "1") == 0;
  }();
  return enabled && neo_event_trace_enabled();
}

CiahCarrierAudit audit_ciah_orbital_geometry(
    const Eigen::VectorXd& orbital_displacement,
    const AcceptedPointContext& origin,
    const StructureAction& action,
    const StructureTangentOperator& structure,
    const std::function<OrbitalPreparationInput(const Eigen::VectorXd&)>& trial_orbitals) {
  CiahCarrierAudit audit;
  if (!action.supports_integral_direction() || !origin.complete_active_space) return audit;
  const auto* inverse = origin.metric_inverse(action);
  if (!inverse) return audit;
  const auto& overlaps = origin.prepared_active_space.orbital_result.active_orbital_overlap_matrix;
  const int n_active = trial_orbitals(Eigen::VectorXd::Zero(orbital_displacement.size()))
      .n_active_orbitals;
  const Eigen::MatrixXd o0 = Eigen::Map<const Eigen::MatrixXd>(overlaps.data(), n_active, n_active);
  const Eigen::MatrixXd r0 = Eigen::LLT<Eigen::MatrixXd>(o0).matrixU();
  const Eigen::MatrixXd raw = action.expand_structure_block(origin.selected_state_eigenvectors);
  const ExteriorOrbitalTransform alpha0(
      origin.same_spin_pair_cache.alpha_reuse_table.unique_determinants, r0);
  const ExteriorOrbitalTransform beta0(
      origin.same_spin_pair_cache.beta_reuse_table.unique_determinants, r0);
  Eigen::MatrixXd y0 = action.orthogonalize_structure_block(origin.selected_state_eigenvectors);
  y0 /= y0.norm();
  struct Geometry {
    Eigen::MatrixXd overlap;
    Eigen::MatrixXd upper;
    Eigen::MatrixXd carrier;
    Eigen::VectorXd gradient;
    double norm = 0.0;
  };
  const auto geometry = [&](double scale) {
    const auto prepared = ActiveSpaceOrbitalPreparer{}.prepare(
        trial_orbitals(scale * orbital_displacement));
    Geometry result;
    result.overlap = Eigen::Map<const Eigen::MatrixXd>(
        prepared.active_orbital_overlap_matrix.data(), n_active, n_active);
    result.upper = Eigen::LLT<Eigen::MatrixXd>(result.overlap).matrixU();
    const ExteriorOrbitalTransform alpha(
        origin.same_spin_pair_cache.alpha_reuse_table.unique_determinants, result.upper);
    const ExteriorOrbitalTransform beta(
        origin.same_spin_pair_cache.beta_reuse_table.unique_determinants, result.upper);
    result.carrier = raw;
    alpha.apply_left(&result.carrier);
    beta.apply_right(&result.carrier);
    result.norm = result.carrier.norm();
    result.carrier /= result.norm;
    // F0 M0^-1 F0^T is identity on this complete selected spin sector.
    // Apply the frozen orthogonal Hamiltonian via the existing origin action;
    // no AO transformation, new CI solve or dense structure matrix is needed.
    Eigen::MatrixXd projected = result.carrier;
    alpha0.apply_adjoint_left(&projected);
    beta0.apply_adjoint_right(&projected);
    const Eigen::MatrixXd coefficients = inverse->apply(
        action.contract_spin_product_block(projected));
    const double closure = (action.orthogonalize_structure_block(coefficients) -
        result.carrier).norm();
    Eigen::MatrixXd sigma = action.orthogonalize_structure_block(inverse->apply(
        action.apply(coefficients).hamiltonian));
    const double energy = (result.carrier.array() * sigma.array()).sum();
    sigma -= energy * result.carrier;
    alpha.apply_adjoint_left(&sigma);
    beta.apply_adjoint_right(&sigma);
    result.gradient = structure.project_coordinates(
        (std::sqrt(2.0 * origin.normalized_state_weights[0]) / result.norm) *
        action.contract_spin_product_block(sigma));
    neo_trace("ciah_carrier_counterfactual", "scale=%.17g,closure=%.17g,gradient=%.17g",
        scale, closure, result.gradient.norm());
    return result;
  };
  const double epsilon = 3.0e-4;
  const Geometry zero = geometry(0.0);
  const Geometry plus = geometry(epsilon), minus = geometry(-epsilon);
  const Eigen::MatrixXd overlap_direction = (plus.overlap - minus.overlap) / (2.0 * epsilon);
  const Eigen::MatrixXd carrier_direction = (plus.carrier - minus.carrier) / (2.0 * epsilon);
  const Eigen::VectorXd gradient_direction = (plus.gradient - minus.gradient) / (2.0 * epsilon);
  audit.remainders.resize(structure.tangent_size(), 4);
  neo_trace("ciah_path_geometry_origin",
      "overlap_identity_error=%.17g,overlap_direction=%.17g,carrier_direction=%.17g",
      (o0 - Eigen::MatrixXd::Identity(n_active, n_active)).norm(),
      overlap_direction.norm(), carrier_direction.norm());
  int column = 0;
  for (const double scale : {0.125, 0.25, 0.5, 1.0}) {
    const Geometry point = geometry(scale);
    audit.remainders.col(column++) = point.gradient - zero.gradient - scale * gradient_direction;
    neo_trace("ciah_path_geometry",
        "scale=%.17g,overlap_change=%.17g,upper_change=%.17g,carrier_change=%.17g,"
        "carrier_nonlinearity=%.17g,overlap_nonlinearity=%.17g,carrier_norm=%.17g",
        scale, (point.overlap - o0).norm(), (point.upper - r0).norm(),
        (point.carrier - y0).norm(),
        (point.carrier - y0 - scale * carrier_direction).norm(),
        (point.overlap - o0 - scale * overlap_direction).norm(), point.norm);
  }
  return audit;
}

void audit_ciah_keyframe(
    const Eigen::VectorXd& displacement,
    const Eigen::VectorXd& origin_gradient,
    const Eigen::VectorXd& hessian_displacement,
    const Eigen::VectorXd& actual_gradient,
    const Eigen::MatrixXd& origin_coefficients,
    double weight,
    const StructureTangentOperator& structure,
    const StructureAction& action,
    const CiahCarrierAudit& carrier,
    const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& evaluate_gradient) {
  const Eigen::Index ns = structure.tangent_size();
  const Eigen::Index no = displacement.size() - ns;
  const Eigen::VectorXd linear = hessian_displacement.tail(ns);
  constexpr double scales[] = {0.125, 0.25, 0.5, 1.0};
  const Eigen::VectorXd delta_ci =
      structure.coefficient_response_coordinates(displacement.tail(ns)).col(0);
  Eigen::MatrixXd ci(origin_coefficients.rows(), 5);
  ci.col(0) = origin_coefficients.col(0);
  for (int j = 0; j < 4; ++j) {
    ci.col(j + 1) = ci.col(0) + scales[j] * delta_ci;
  }
  // Measure the accepted CI residual rather than assuming it is exactly zero.
  const StructureActionResult images = action.apply(ci);
  const auto ci_gradient = [&](int column) -> Eigen::VectorXd {
    const double norm_squared = ci.col(column).dot(images.overlap.col(column));
    const double energy = ci.col(column).dot(images.hamiltonian.col(column)) / norm_squared;
    const Eigen::MatrixXd residual = (std::sqrt(2.0 * weight) / norm_squared) *
        (images.hamiltonian.col(column) - energy * images.overlap.col(column));
    return structure.project_coordinates(residual);
  };
  const Eigen::VectorXd g0 = ci_gradient(0);
  neo_trace("ciah_path_origin",
      "structure_gradient=%.17g,structure_gradient_physical=%.17g,forced_zero_error=%.17g",
      g0.norm(), structure.dual_norm(g0), (g0 - origin_gradient.tail(ns)).norm());

  // Test the actual accumulated path, not another arbitrary endpoint probe.
  for (const double epsilon : {1.0e-3, 3.0e-4}) {
    const Eigen::VectorXd plus = evaluate_gradient(epsilon * displacement);
    const Eigen::VectorXd minus = evaluate_gradient(-epsilon * displacement);
    const Eigen::VectorXd finite = (plus - minus) / (2.0 * epsilon);
    const Eigen::VectorXd error = finite - hessian_displacement;
    neo_trace("ciah_path_derivative",
        "epsilon=%.17g,orbital_error=%.17g,structure_error=%.17g,"
        "orbital_image=%.17g,structure_image=%.17g",
        epsilon, error.head(no).norm(), error.tail(ns).norm(),
        hessian_displacement.head(no).norm(), linear.norm());
  }

  Eigen::VectorXd orbital_displacement = displacement;
  orbital_displacement.tail(ns).setZero();
  const Eigen::VectorXd ci_linear =
      structure.apply_hessian_coordinates(displacement.tail(ns));
  const Eigen::VectorXd orbital_linear = linear - ci_linear;
  const double epsilon = 3.0e-4;
  const Eigen::VectorXd orbital_plus = evaluate_gradient(epsilon * orbital_displacement);
  const Eigen::VectorXd orbital_minus = evaluate_gradient(-epsilon * orbital_displacement);
  const Eigen::VectorXd orbital_finite =
      (orbital_plus.tail(ns) - orbital_minus.tail(ns)) / (2.0 * epsilon);
  neo_trace("ciah_path_coupling_derivative",
      "epsilon=%.17g,structure_error=%.17g,structure_image=%.17g",
      epsilon, (orbital_finite - orbital_linear).norm(), orbital_linear.norm());
  for (int j = 0; j < 4; ++j) {
    const double scale = scales[j];
    const Eigen::VectorXd orbital_gradient =
        evaluate_gradient(scale * orbital_displacement).tail(ns);
    const Eigen::VectorXd ci_g = ci_gradient(j + 1);
    const Eigen::VectorXd actual = scale == 1.0 ? Eigen::VectorXd(actual_gradient.tail(ns)) :
        Eigen::VectorXd(evaluate_gradient(scale * displacement).tail(ns));
    const Eigen::VectorXd orbital_error = orbital_gradient - g0 - scale * orbital_linear;
    const Eigen::VectorXd ci_error = ci_g - g0 - scale * ci_linear;
    const Eigen::VectorXd mixed = actual - orbital_gradient - ci_g + g0;
    const Eigen::VectorXd error = actual - g0 - scale * linear;
    if (carrier.remainders.cols() == 4) {
      const auto representation = carrier.remainders.col(j);
      neo_trace("ciah_path_orbital_split",
          "scale=%.17g,representation_remainder=%.17g,remaining_remainder=%.17g,"
          "total_remainder=%.17g,cross=%.17g",
          scale, representation.norm(), (orbital_error - representation).norm(),
          orbital_error.norm(), representation.dot(orbital_error - representation));
    }
    const double norm_squared = ci.col(j + 1).dot(images.overlap.col(j + 1));
    neo_trace("ciah_path_structure",
      "scale=%.17g,"
      "orbital_linear=%.17g,ci_linear=%.17g,orbital_actual=%.17g,ci_actual=%.17g,"
      "orbital_error=%.17g,ci_error=%.17g,mixed=%.17g,total_error=%.17g,"
      "orbital_error_physical=%.17g,ci_error_physical=%.17g,mixed_physical=%.17g,"
      "total_error_physical=%.17g,reconstruction=%.17g,ci_norm_squared=%.17g,"
      "orbital_ci_error_dot=%.17g,orbital_mixed_dot=%.17g,ci_mixed_dot=%.17g",
      scale, scale * orbital_linear.norm(), scale * ci_linear.norm(),
      orbital_gradient.norm(), ci_g.norm(),
      orbital_error.norm(), ci_error.norm(), mixed.norm(), error.norm(),
      structure.dual_norm(orbital_error), structure.dual_norm(ci_error),
      structure.dual_norm(mixed), structure.dual_norm(error),
      (orbital_error + ci_error + mixed - error).norm(), norm_squared,
      orbital_error.dot(ci_error), orbital_error.dot(mixed), ci_error.dot(mixed));
  }
}

}  // namespace xmvb::vb::diagnostics
