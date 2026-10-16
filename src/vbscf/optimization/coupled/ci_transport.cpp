#include "vbscf/optimization/coupled/ci_transport.hpp"

#include <cmath>
#include <chrono>
#include "core/metric_solve.hpp"
#include "core/timing.hpp"
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
      c0_(origin.selected_state_eigenvectors),
      coordinate_scale_(std::sqrt(2.0 * origin.normalized_state_weights.at(0))),
      diagonal_(action.preconditioner_diagonal().overlap),
      inverse_(origin.metric_inverse(action)) {
  const auto& images = origin.selected_state_structure_images.value();
  const Eigen::MatrixXd state_h = c0_.transpose() * images.hamiltonian;
  gc0_ = (coordinate_scale_ * coordinate_scale_) *
      (images.hamiltonian - images.overlap * state_h);
  lambda0_ = solve(upper_, gc0_);
  residual_upper_gradient_ = zero_;
  for (Eigen::Index s = 0; s < c0_.cols(); ++s)
    residual_upper_gradient_ += overlap_.adjoint(upper_, upper_, lambda0_.col(s),
        c0_.col(s),OverlapSide::Right).right;
  residual_overlap_gradient_ = overlap_upper_gradient(upper_, residual_upper_gradient_);
}

Eigen::VectorXd CiTransport::raw_direction(const Eigen::VectorXd& reduced) const {
  return gather_nonredundant_retract_tangent(input_, chart_, layout_, reduced);
}

Eigen::MatrixXd CiTransport::solve(const Eigen::MatrixXd& upper,
                                 const Eigen::MatrixXd& rhs,
                                 const Eigen::MatrixXd* initial) const {
  // Opt-in diagnosis only: separate the selected metric solve from HVP work.
  const bool timed = profiling::phase_totals().enabled;
  using Clock = std::chrono::steady_clock;
  const auto start = timed ? Clock::now() : Clock::time_point{};
  double action_seconds = 0.0;
  const auto report = [&](const char* mode, int actions, int iterations,
                          int certifications, double residual) {
    if (timed)
      std::fprintf(stderr,
          "XMVB_TRANSPORT_SOLVE,mode=%s,rows=%ld,rhs=%ld,actions=%d,iterations=%d,certifications=%d,relative_residual=%.9g,wall=%.9f,action_wall=%.9f\n",
          mode, static_cast<long>(rhs.rows()), static_cast<long>(rhs.cols()),
          actions, iterations, certifications, residual,
          std::chrono::duration<double>(Clock::now() - start).count(), action_seconds);
  };
  if (inverse_) {
    const ExteriorOrbitalTransform alpha(
        origin_.same_spin_pair_cache.alpha_reuse_table.unique_determinants, upper);
    const ExteriorOrbitalTransform beta(
        origin_.same_spin_pair_cache.beta_reuse_table.unique_determinants, upper);
    Eigen::MatrixXd result = inverse_->apply_at_transform(rhs, alpha, beta);
    report("factorized", 0, 0, 0, 0.0);
    return result;
  }
  // Actual selected M, not the identity on a presumed closed carrier.
  if (!prepared_metric_ || metric_upper_.rows() != upper.rows() ||
      !(metric_upper_.array() == upper.array()).all()) {
    prepared_metric_.emplace(overlap_.prepare_metric(upper));
    metric_upper_ = upper;
  }
  // Solve only the selected Gram. The positive accepted-point diagonal is a
  // Jacobi preconditioner also at a trial chart: it changes convergence, not
  // the normal equation. Do not complete the selected space to precondition it.
  Eigen::MatrixXd result(rhs.rows(), rhs.cols());
  int actions = 0, iterations = 0, certifications = 0;
  double residual = 0.0;
  for (Eigen::Index s = 0; s < rhs.cols(); ++s) {
    core::MetricSolveStats stats;
    Eigen::VectorXd seed;
    if (initial) seed = initial->col(s);
    result.col(s) = core::solve_metric(rhs.col(s),
      [&](const Eigen::VectorXd& x) {
        const auto action_start = timed ? Clock::now() : Clock::time_point{};
        Eigen::VectorXd image = (*prepared_metric_)(x).col(0);
        if (timed) action_seconds += std::chrono::duration<double>(Clock::now() - action_start).count();
        return image;
      },
      [&](const Eigen::VectorXd& x) {
        return Eigen::VectorXd(x.cwiseQuotient(diagonal_));
      },
      timed ? &stats : nullptr, initial ? &seed : nullptr);
    actions += stats.actions;
    iterations += stats.iterations;
    certifications += stats.certifications;
    residual = std::max(residual, stats.relative_residual);
  }
  report("selected_cg", actions, iterations, certifications, residual);
  return result;
}

CiTransportPoint CiTransport::retract(const Eigen::VectorXd& orbital_step,
                                    const Eigen::VectorXd& structure_step) const {
  profiling::ScopedPhase phase(profiling::Phase::CiTransportRetract);
  CiTransportPoint point;
  point.input = chart_.retract_step(input_, orbital_step);
  point.prepared = ActiveSpaceOrbitalPreparer{}.prepare(point.input);
  point.upper = overlap_upper(Eigen::Map<const Eigen::MatrixXd>(
      point.prepared.active_orbital_overlap_matrix.data(), input_.n_active_orbitals,
      input_.n_active_orbitals));
  point.origin_direction = overlap_upper_direction(upper_,
      origin_pullback_.overlap_direction(raw_direction(orbital_step)));
  point.base_coefficients = c0_ +
      structure_.coefficient_response_coordinates(structure_step);
  const Eigen::MatrixXd rhs = overlap_.apply(point.upper, upper_, point.base_coefficients) +
      overlap_.apply_direction(point.upper, upper_, zero_, point.origin_direction, c0_);
  // The target already agrees with F' * base_coefficients to first order
  // in the orbital step (at zero independent CI step). Correct the defect
  // of this known representative, with the original RHS-relative tolerance.
  point.coefficients = solve(point.upper, rhs, &point.base_coefficients);
  return point;
}

Eigen::VectorXd CiTransport::gradient(const CiTransportPoint& point,
                                    const VbScfCoupledKeyframeEvaluation& fixed) const {
  profiling::ScopedPhase phase(profiling::Phase::CiTransportGradient);
  const Eigen::MatrixXd lambda = solve(point.upper, fixed.coefficient_gradient);
  Eigen::MatrixXd current_upper = zero_, origin_upper = zero_;
  Eigen::MatrixXd structure_gradient(c0_.rows(),c0_.cols());
  for (Eigen::Index s = 0; s < c0_.cols(); ++s) {
    const auto pullback=overlap_.transport_gradient(point.upper,upper_,point.origin_direction,
        lambda.col(s),point.base_coefficients.col(s),c0_.col(s),point.coefficients.col(s));
    current_upper+=pullback.charts.left;
    origin_upper+=pullback.charts.right;
    structure_gradient.col(s)=pullback.image/coordinate_scale_;
  }
  const OrbitalOverlapPullback current_pullback(point.input, point.prepared, layout_);
  Eigen::VectorXd orbital = fixed.orbital_gradient +
      current_pullback.gradient(overlap_upper_gradient(point.upper, current_upper)) +
      origin_pullback_.gradient(overlap_upper_gradient(upper_, origin_upper));
  Eigen::VectorXd result(chart_.reduced_size() + structure_.tangent_size());
  result.head(chart_.reduced_size()) = chart_.project_reduced_gradient(orbital);
  result.tail(structure_.tangent_size()) = structure_.project_coordinates(structure_gradient);
  return result;
}

CiTransportCurvature CiTransport::curvature(const Eigen::VectorXd& orbital_direction,
                                          const Eigen::VectorXd& structure_direction) const {
  profiling::ScopedPhase phase(profiling::Phase::CiTransportCurvature);
  const Eigen::VectorXd packed = raw_direction(orbital_direction);
  const Eigen::MatrixXd dr = overlap_upper_direction(upper_,
      origin_pullback_.overlap_direction(packed));
  const Eigen::MatrixXd dc = structure_.coefficient_response_coordinates(structure_direction);
  Eigen::MatrixXd dg = zero_, mixed = zero_;
  Eigen::MatrixXd forcing(c0_.rows(),c0_.cols());
  for (Eigen::Index s = 0; s < c0_.cols(); ++s) {
    const auto contractions=overlap_.transport_curvature(upper_,dr,lambda0_.col(s),
        c0_.col(s),dc.col(s));
    dg+=contractions.directional;
    mixed+=contractions.mixed;
    forcing.col(s)=contractions.image;
  }
  const Eigen::MatrixXd d_overlap_gradient = overlap_upper_gradient_direction(
      upper_, dr, residual_upper_gradient_, dg);
  CiTransportCurvature result;
  result.orbital_from_orbital = -chart_.project_reduced_gradient(
      origin_pullback_.gradient_direction(residual_overlap_gradient_, d_overlap_gradient, packed));
  result.orbital_from_structure = -chart_.project_reduced_gradient(
      origin_pullback_.gradient(overlap_upper_gradient(upper_, mixed)));
  result.structure_from_orbital = -structure_.project_coordinates(
      forcing / coordinate_scale_);
  return result;
}

Eigen::VectorXd CiTransport::origin_structure_gradient() const {
  return structure_.project_coordinates(gc0_ / coordinate_scale_);
}

}  // namespace xmvb::vb
