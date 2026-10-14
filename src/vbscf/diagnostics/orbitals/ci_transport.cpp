#include "vbscf/diagnostics/orbitals/ci_transport.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <Eigen/Cholesky>

#include "vbscf/derivatives/hessian/context/accepted_point.hpp"
#include "vbscf/derivatives/hessian/coupled/structure.hpp"
#include "vbscf/derivatives/hessian/responses/orbital/preparation.hpp"
#include "vbscf/diagnostics/performance/neo_trace.hpp"
#include "vbscf/optimization/objective/function.hpp"
#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"
#include "vbscf/structures/orthogonal_ci/metric.hpp"

namespace xmvb::vb::diagnostics {
namespace {

struct Geometry {
  std::vector<double> overlap;
  Eigen::MatrixXd upper;
  Eigen::MatrixXd direction;
};

/** @brief Overlap-only derivative along the immutable raw orbital tangent. */
Geometry geometry(const OrbitalPreparationInput& input,
                  const SparseParameterLayout& layout,
                  const Eigen::VectorXd& raw_direction) {
  const auto prepared = ActiveSpaceOrbitalPreparer{}.prepare(input);
  const auto cache = build_accepted_orbital_preparation_cache(input, Eigen::MatrixXd{}, {});
  const auto tangent = build_dense_orbital_tangent_context(input, layout, raw_direction, cache);
  const int nc = (input.n_total_electrons - input.n_active_electrons) / 2;
  const int na = input.n_active_orbitals;
  const auto direction = build_orbital_preparation_directional_result(
      input, tangent, nc, na, cache);
  const Eigen::Map<const Eigen::MatrixXd> overlap(
      input.ao_overlap_matrix.data(), input.n_basis_functions, input.n_basis_functions);
  const Eigen::MatrixXd left = direction.delta_active_auxiliary_orbitals.transpose() *
      overlap * prepared.auxiliary_orbital_matrix.middleCols(nc, na);
  const Eigen::MatrixXd delta = left + left.transpose();
  Geometry result;
  result.overlap = prepared.active_orbital_overlap_matrix;
  const Eigen::Map<const Eigen::MatrixXd> active_overlap(result.overlap.data(), na, na);
  Eigen::LLT<Eigen::MatrixXd> factor(active_overlap);
  if (factor.info() != Eigen::Success) throw std::runtime_error("transport overlap is singular");
  result.upper = factor.matrixU();
  const Eigen::MatrixXd inverse = result.upper.triangularView<Eigen::Upper>().solve(
      Eigen::MatrixXd::Identity(na, na));
  const Eigen::MatrixXd normalized = inverse.transpose() * delta * inverse;
  Eigen::MatrixXd generator = normalized.triangularView<Eigen::Upper>();
  generator.diagonal() *= 0.5;
  result.direction = generator * result.upper;
  return result;
}

/** @brief Carrier maps share the origin structure expansion; no closure assumption. */
class CarrierMap {
public:
  CarrierMap(const Geometry& g, const AcceptedPointContext& origin, const StructureAction& action)
      : action_(action),
        alpha_(origin.same_spin_pair_cache.alpha_reuse_table.unique_determinants, g.upper),
        beta_(origin.same_spin_pair_cache.beta_reuse_table.unique_determinants, g.upper),
        alpha_direction_(alpha_.direction(g.direction)),
        beta_direction_(beta_.direction(g.direction)) {}

  Eigen::MatrixXd apply(const Eigen::MatrixXd& c) const {
    Eigen::MatrixXd y = action_.expand_structure_block(c);
    alpha_.apply_left(&y);
    beta_.apply_right(&y);
    return y;
  }

  Eigen::MatrixXd adjoint(Eigen::MatrixXd y) const {
    beta_.apply_adjoint_right(&y);
    alpha_.apply_adjoint_left(&y);
    return action_.contract_spin_product_block(y);
  }

  Eigen::MatrixXd direction(const Eigen::MatrixXd& c) const {
    Eigen::MatrixXd y = action_.expand_structure_block(c);
    Eigen::MatrixXd dy = Eigen::MatrixXd::Zero(y.rows(), y.cols());
    alpha_.apply_directional_left(alpha_direction_, &y, &dy);
    beta_.apply_directional_right(beta_direction_, &y, &dy);
    return dy;
  }

  Eigen::MatrixXd direction_adjoint(Eigen::MatrixXd y) const {
    Eigen::MatrixXd dy = Eigen::MatrixXd::Zero(y.rows(), y.cols());
    beta_.apply_directional_adjoint_right(beta_direction_, &y, &dy);
    alpha_.apply_directional_adjoint_left(alpha_direction_, &y, &dy);
    return action_.contract_spin_product_block(dy);
  }

private:
  const StructureAction& action_;
  ExteriorOrbitalTransform alpha_, beta_;
  ExteriorTransformDirection alpha_direction_, beta_direction_;
};

struct Point {
  std::unique_ptr<StructureAction> action;
  std::unique_ptr<StructureMetricInverse> inverse;
  Eigen::MatrixXd coefficients;
  Eigen::MatrixXd residual;
  Eigen::MatrixXd coefficient_direction;
  Eigen::VectorXd structure_gradient;
  double energy = 0.0;
  double joint_derivative = 0.0;
};

}  // namespace

bool ciah_ci_transport_audit_enabled() {
  const char* flag = std::getenv("XMVB_NEO_KEYFRAME_AUDIT");
  return neo_event_trace_enabled() && flag && std::strcmp(flag, "transport") == 0;
}

void audit_ciah_ci_transport(
    const Eigen::VectorXd& displacement,
    const Eigen::VectorXd& predicted_structure,
    const Eigen::VectorXd& actual_structure,
    const AcceptedPointContext& origin,
    const FullDeterminantStructureData& structures,
    const StructureAction& action,
    const StructureTangentOperator& structure,
    const SparseParameterLayout& layout,
    const std::function<OrbitalPreparationInput(const Eigen::VectorXd&)>& trial_orbitals,
    const std::function<VbScfCoupledKeyframeEvaluation(
        const Eigen::VectorXd&, const Eigen::MatrixXd&)>& evaluate) {
  // The projection theory is general. This molecular prototype deliberately
  // tests a complete sector; nonclosed selected sectors have a separate oracle.
  if (!origin.complete_active_space || !action.supports_integral_direction())
    throw std::invalid_argument("molecular transport audit requires the admitted complete sector");
  const auto start = std::chrono::steady_clock::now();
  const Eigen::Index ns = structure.tangent_size(), no = displacement.size() - ns;
  const Eigen::VectorXd p = displacement.head(no), q = displacement.tail(ns);
  const OrbitalPreparationInput input = trial_orbitals(Eigen::VectorXd::Zero(no));
  const Eigen::VectorXd raw_direction = layout.pack(trial_orbitals(p)) - layout.pack(input);
  const Geometry g0 = geometry(input, layout, raw_direction);
  const CarrierMap f0(g0, origin, action);
  const Eigen::MatrixXd c0 = origin.selected_state_eigenvectors;
  const Eigen::MatrixXd dq = structure.coefficient_response_coordinates(q);
  const Eigen::MatrixXd connection = f0.direction(c0);
  const auto* inverse0 = origin.metric_inverse(action);
  if (!inverse0) throw std::runtime_error("origin transport metric factors are unavailable");

  const auto point = [&](double orbital_scale, double ci_scale, bool gradient) {
    const auto trial_input = trial_orbitals(orbital_scale * p);
    const Geometry g = geometry(trial_input, layout, raw_direction);
    const CarrierMap f(g, origin, action);
    Point result;
    const auto& prepared = origin.prepared_active_space;
    result.action = std::make_unique<StructureAction>(
        structures.determinant_to_structure_terms, action.n_structures(), origin.same_spin_pair_cache,
        g.overlap, prepared.active_space_one_electron_result.h1e_act,
        prepared.active_space_two_electron_result, input.n_active_orbitals,
        nullptr, StructureActionMode::ActionOnly);
    // Own trial factors; never rebind the immutable origin inverse.
    result.inverse = std::make_unique<StructureMetricInverse>(*result.action,
        origin.same_spin_pair_cache.alpha_reuse_table.unique_determinants,
        origin.same_spin_pair_cache.beta_reuse_table.unique_determinants, true,
        64 * 1024 * 1024, 64 * 1024 * 1024);
    const Eigen::MatrixXd target = f0.apply(c0 + ci_scale * dq) + orbital_scale * connection;
    const Eigen::MatrixXd rhs = f.adjoint(target);
    result.coefficients = result.inverse->apply(rhs);
    result.residual = target - f.apply(result.coefficients);
    const Eigen::MatrixXd target_direction = f0.apply(dq) + connection;
    const Eigen::MatrixXd dc_rhs = f.direction_adjoint(result.residual) +
        f.adjoint(target_direction - f.direction(result.coefficients));
    result.coefficient_direction = result.inverse->apply(dc_rhs);
    if (gradient) {
      const auto exact = evaluate(orbital_scale * p, result.coefficients);
      result.energy = exact.energy;
      const double norm = std::sqrt(result.coefficients.col(0).dot(
          result.action->apply_overlap(result.coefficients).col(0)));
      const Eigen::MatrixXd gc = (2.0 / norm) * exact.structure_residuals;
      const Eigen::MatrixXd lambda = result.inverse->apply(gc);
      result.structure_gradient = structure.project_coordinates(
          f0.adjoint(f.apply(lambda)) / std::sqrt(2.0));
      result.joint_derivative = exact.orbital_gradient.dot(raw_direction) +
          gc.col(0).dot(result.coefficient_direction.col(0));
    }
    neo_trace("ciah_transport_point",
        "orbital_scale=%.17g,ci_scale=%.17g,projection_residual=%.17g,"
        "normal_residual=%.17g,metric_residual=%.17g,coefficient_correction=%.17g",
        orbital_scale, ci_scale, result.residual.norm(), f.adjoint(result.residual).norm(),
        (result.action->apply_overlap(result.coefficients) - rhs).norm() / rhs.norm(),
        (result.coefficients - c0 - ci_scale * dq).norm());
    return result;
  };

  const Point finite = point(1.0, 1.0, true);
  neo_trace("ciah_transport_response",
      "old_error=%.17g,new_error=%.17g,old_gradient=%.17g,new_gradient=%.17g,"
      "predicted_gradient=%.17g,energy=%.17g,joint_derivative=%.17g",
      (actual_structure - predicted_structure).norm(),
      (finite.structure_gradient - predicted_structure).norm(), actual_structure.norm(),
      finite.structure_gradient.norm(), predicted_structure.norm(), finite.energy,
      finite.joint_derivative);

  // Exact CI-residual correction for the q block of the origin joint HVP.
  const auto& selected = origin.selected_state_structure_images.value();
  const Eigen::MatrixXd gc0 = 2.0 * (selected.hamiltonian -
      origin.selected_state_energies[0] * selected.overlap);
  const Eigen::VectorXd correction = -structure.project_coordinates(
      f0.direction_adjoint(f0.apply(inverse0->apply(gc0))) / std::sqrt(2.0));
  for (const double epsilon : {1.0e-3, 3.0e-4}) {
    const Point plus = point(epsilon, epsilon, true);
    const Point minus = point(-epsilon, -epsilon, true);
    const Eigen::VectorXd image = (plus.structure_gradient - minus.structure_gradient) /
        (2.0 * epsilon);
    neo_trace("ciah_transport_hvp",
        "epsilon=%.17g,structure_error=%.17g,residual_correction=%.17g,"
        "first_jet_error=%.17g",
        epsilon, (image - predicted_structure - correction).norm(), correction.norm(),
        ((plus.coefficients - minus.coefficients) / (2.0 * epsilon) - dq).norm());
  }
  const double epsilon = 3.0e-4;
  const Point plus = point(1.0 + epsilon, 1.0 + epsilon, true);
  const Point minus = point(1.0 - epsilon, 1.0 - epsilon, true);
  neo_trace("ciah_transport_joint_derivative",
      "analytic=%.17g,finite=%.17g,error=%.17g,coefficient_error=%.17g",
      finite.joint_derivative, (plus.energy - minus.energy) / (2.0 * epsilon),
      std::abs(finite.joint_derivative - (plus.energy - minus.energy) / (2.0 * epsilon)),
      ((plus.coefficients - minus.coefficients) / (2.0 * epsilon) -
          finite.coefficient_direction).norm());
  const Point ci_plus = point(1.0, 1.0 + epsilon, true);
  const Point ci_minus = point(1.0, 1.0 - epsilon, true);
  neo_trace("ciah_transport_ci_derivative", "analytic=%.17g,finite=%.17g,error=%.17g",
      finite.structure_gradient.dot(q), (ci_plus.energy - ci_minus.energy) / (2.0 * epsilon),
      std::abs(finite.structure_gradient.dot(q) -
          (ci_plus.energy - ci_minus.energy) / (2.0 * epsilon)));
  neo_trace("ciah_transport_audit_cost", "seconds=%.17g",
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
}

}  // namespace xmvb::vb::diagnostics
