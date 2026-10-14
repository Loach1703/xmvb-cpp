#include "vbscf/optimization/coupled/ci_transport.hpp"

#include <cmath>
#include "core/metric_solve.hpp"
#include "vbscf/derivatives/hessian/context/accepted_point.hpp"
#include "vbscf/derivatives/hessian/coupled/structure.hpp"
#include "vbscf/diagnostics/performance/neo_trace.hpp"
#include "vbscf/optimization/objective/function.hpp"
#include "vbscf/optimization/trust_region/retraction.hpp"
#include "vbscf/orbitals/preparation/cholesky.hpp"
#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"
#include "vbscf/structures/orthogonal_ci/metric.hpp"

namespace xmvb::vb {

CiTransport::CiTransport(const OrbitalPreparationInput& input,
                         const SparseParameterLayout& layout,
                         const OrbitalChart& chart,
                         const AcceptedPointContext& origin,
                         const StructureAction& action,
                         const StructureTangentOperator& structure)
    : input_(input), layout_(layout), chart_(chart), origin_(origin), action_(action),
      structure_(structure),
      overlap_(action, origin.same_spin_pair_cache.alpha_reuse_table.unique_determinants,
          origin.same_spin_pair_cache.beta_reuse_table.unique_determinants, input.n_active_orbitals),
      origin_pullback_(input, origin.prepared_active_space.orbital_result, layout),
      upper_(overlap_upper(Eigen::Map<const Eigen::MatrixXd>(
          origin.prepared_active_space.orbital_result.active_orbital_overlap_matrix.data(),
          input.n_active_orbitals, input.n_active_orbitals))),
      zero_(Eigen::MatrixXd::Zero(input.n_active_orbitals, input.n_active_orbitals)),
      c0_(origin.selected_state_eigenvectors.col(0)),
      diagonal_(action.preconditioner_diagonal().overlap),
      inverse_(origin.metric_inverse(action)) {
  if (origin.selected_state_indices.size() != 1)
    throw std::invalid_argument("finite CI transport requires one selected state");
  const auto& images = origin.selected_state_structure_images.value();
  gc0_ = 2.0 * (images.hamiltonian.col(0) -
      origin.selected_state_energies[0] * images.overlap.col(0));
  lambda0_ = solve(upper_, gc0_);
  residual_upper_gradient_ = overlap_.adjoint(upper_, upper_, lambda0_, c0_).right;
  residual_overlap_gradient_ = overlap_upper_gradient(upper_, residual_upper_gradient_);
}

Eigen::VectorXd CiTransport::raw_direction(const Eigen::VectorXd& reduced) const {
  return gather_nonredundant_retract_tangent(input_, chart_, layout_, reduced);
}

Eigen::VectorXd CiTransport::solve(const Eigen::MatrixXd& upper,
                                 const Eigen::VectorXd& rhs) const {
  if (inverse_) {
    const ExteriorOrbitalTransform alpha(
        origin_.same_spin_pair_cache.alpha_reuse_table.unique_determinants, upper);
    const ExteriorOrbitalTransform beta(
        origin_.same_spin_pair_cache.beta_reuse_table.unique_determinants, upper);
    return inverse_->apply_at_transform(rhs, alpha, beta).col(0);
  }
  // Actual selected M, not the identity on a presumed closed carrier.
  return core::solve_metric(rhs,
      [&](const Eigen::VectorXd& x) {
        return Eigen::VectorXd(overlap_.apply(upper, upper, x).col(0));
      },
      [&](const Eigen::VectorXd& x) { return Eigen::VectorXd(x.cwiseQuotient(diagonal_)); });
}

CiTransportPoint CiTransport::retract(const Eigen::VectorXd& orbital_step,
                                    const Eigen::VectorXd& structure_step) const {
  CiTransportPoint point;
  point.input = chart_.retract_step(input_, orbital_step);
  point.prepared = ActiveSpaceOrbitalPreparer{}.prepare(point.input);
  point.upper = overlap_upper(Eigen::Map<const Eigen::MatrixXd>(
      point.prepared.active_orbital_overlap_matrix.data(), input_.n_active_orbitals,
      input_.n_active_orbitals));
  point.origin_direction = overlap_upper_direction(upper_,
      origin_pullback_.overlap_direction(raw_direction(orbital_step)));
  point.base_coefficients = c0_ +
      structure_.coefficient_response_coordinates(structure_step).col(0);
  const Eigen::VectorXd rhs = overlap_.apply(point.upper, upper_, point.base_coefficients).col(0) +
      overlap_.apply_direction(point.upper, upper_, zero_, point.origin_direction, c0_).col(0);
  point.coefficients = solve(point.upper, rhs);
  return point;
}

Eigen::VectorXd CiTransport::gradient(const CiTransportPoint& point,
                                    const VbScfCoupledKeyframeEvaluation& fixed) const {
  const Eigen::VectorXd c = point.coefficients.col(0);
  const double norm = std::sqrt(c.dot(overlap_.apply(point.upper, point.upper, c).col(0)));
  const Eigen::VectorXd gc = (2.0 / norm) * fixed.structure_residuals.col(0);
  const Eigen::VectorXd lambda = solve(point.upper, gc);
  const UpperPair base = overlap_.adjoint(point.upper, upper_, lambda, point.base_coefficients);
  const UpperPair origin = overlap_.adjoint(point.upper, upper_, lambda, c0_);
  const UpperPair connection = overlap_.adjoint_direction(point.upper, upper_,
      zero_, point.origin_direction, lambda, c0_);
  const UpperPair metric = overlap_.adjoint(point.upper, point.upper, lambda, c);
  const Eigen::MatrixXd current_upper = base.left + connection.left - metric.left - metric.right;
  const OrbitalOverlapPullback current_pullback(point.input, point.prepared, layout_);
  Eigen::VectorXd orbital = fixed.orbital_gradient +
      current_pullback.gradient(overlap_upper_gradient(point.upper, current_upper)) +
      origin_pullback_.gradient(overlap_upper_gradient(upper_, origin.right));
  Eigen::VectorXd result(chart_.reduced_size() + structure_.tangent_size());
  result.head(chart_.reduced_size()) = chart_.project_reduced_gradient(orbital);
  result.tail(structure_.tangent_size()) = structure_.project_coordinates(
      overlap_.apply(upper_, point.upper, lambda) / std::sqrt(2.0));
  return result;
}

CiTransportCurvature CiTransport::curvature(const Eigen::VectorXd& orbital_direction,
                                          const Eigen::VectorXd& structure_direction) const {
  const Eigen::VectorXd packed = raw_direction(orbital_direction);
  const Eigen::MatrixXd dr = overlap_upper_direction(upper_,
      origin_pullback_.overlap_direction(packed));
  const Eigen::VectorXd dc = structure_.coefficient_response_coordinates(structure_direction).col(0);
  const Eigen::MatrixXd dg = overlap_.adjoint_direction(
      upper_, upper_, zero_, dr, lambda0_, c0_).right;
  const Eigen::MatrixXd d_overlap_gradient = overlap_upper_gradient_direction(
      upper_, dr, residual_upper_gradient_, dg);
  CiTransportCurvature result;
  result.orbital_from_orbital = -chart_.project_reduced_gradient(
      origin_pullback_.gradient_direction(residual_overlap_gradient_, d_overlap_gradient, packed));
  const Eigen::MatrixXd mixed = overlap_.adjoint(upper_, upper_, lambda0_, dc).right;
  result.orbital_from_structure = -chart_.project_reduced_gradient(
      origin_pullback_.gradient(overlap_upper_gradient(upper_, mixed)));
  result.structure_from_orbital = -structure_.project_coordinates(
      overlap_.apply_direction(upper_, upper_, dr, zero_, lambda0_) / std::sqrt(2.0));
  return result;
}

Eigen::VectorXd CiTransport::origin_structure_gradient() const {
  return structure_.project_coordinates(Eigen::MatrixXd(gc0_ / std::sqrt(2.0)));
}

}  // namespace xmvb::vb
