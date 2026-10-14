#include "vbscf/diagnostics/orbitals/step_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

#include "vbscf/diagnostics/performance/neo_trace.hpp"
#include "vbscf/orbitals/charts/chart.hpp"
#include "vbscf/orbitals/preparation/input.hpp"

namespace xmvb::vb::diagnostics {

OrbitalStepGeometry measure_orbital_step(
    const Eigen::Ref<const Eigen::VectorXd>& coefficients,
    const Eigen::Ref<const Eigen::VectorXd>& displacement,
    const Eigen::Ref<const Eigen::MatrixXd>& overlap) {
  if (coefficients.size() != displacement.size() ||
      overlap.rows() != coefficients.size() || overlap.cols() != coefficients.size())
    throw std::invalid_argument("orbital step geometry dimensions disagree");
  const double norm = std::sqrt(coefficients.dot(overlap * coefficients));
  if (!(norm > 0.0) || !std::isfinite(norm))
    throw std::invalid_argument("orbital step geometry has invalid origin norm");
  const Eigen::VectorXd x = coefficients / norm;
  const Eigen::VectorXd d = displacement / norm;
  OrbitalStepGeometry result;
  result.alpha = x.dot(overlap * d);
  const Eigen::VectorXd tangent = d - result.alpha * x;
  // Form the tangent norm directly: ||d||^2-alpha^2 loses accuracy for
  // nearly radial displacements. No AO matrix is copied or factored here.
  result.theta = std::sqrt(std::max(0.0, tangent.dot(overlap * tangent)));
  result.step_norm = std::hypot(result.alpha, result.theta);
  result.norm_ratio = std::hypot(1.0 + result.alpha, result.theta);
  result.actual_angle = std::atan2(result.theta, 1.0 + result.alpha);
  return result;
}

bool orbital_step_audit_enabled() {
  static const bool enabled = [] {
    const char* flag = std::getenv("XMVB_NEO_ORBITAL_STEP_AUDIT");
    return flag && std::strcmp(flag, "1") == 0;
  }();
  return enabled && neo_event_trace_enabled();
}

void audit_orbital_step(const OrbitalPreparationInput& input,
                        const OrbitalChart& chart,
                        const Eigen::VectorXd& reduced_step,
                        const char* stage, int micro,
                        const Eigen::VectorXd& base_reduced) {
  const Eigen::VectorXd raw = chart.expand_retract_input_tangent(input, reduced_step);
  Eigen::VectorXd offset;
  if (base_reduced.size()) offset = chart.expand_retract_input_tangent(input, base_reduced);
  const int n_inactive = (input.n_total_electrons - input.n_active_electrons) / 2;
  const int n_occupied = n_inactive + input.n_active_orbitals;
  // Two AO vectors are reused serially. The shared AO overlap is borrowed;
  // neither a per-orbital metric copy nor a per-thread workspace is needed.
  Eigen::VectorXd b(input.n_basis_functions), delta(input.n_basis_functions);
  for (int orbital = 0; orbital < n_occupied; ++orbital) {
    b.setZero();
    delta.setZero();
    const int count = stored_sparse_orbital_coefficient_count(input, orbital);
    for (int slot = 0; slot < count; ++slot) {
      const auto flat = orbital * input.n_basis_functions + slot;
      const int basis = input.orbital_basis_index_table[flat] - 1;
      b[basis] = input.orbital_value_table[flat] + (offset.size() ? offset[flat] : 0.0);
      delta[basis] = raw[flat];
    }
    const auto geometry = measure_orbital_step(b, delta, input.ao_overlap_matrix);
    neo_trace("orbital_step_geometry",
        "stage=%s,micro=%d,orbital=%d,role=%s,coefficients=%d,alpha=%.17g,theta=%.17g,"
        "step_norm=%.17g,actual_angle=%.17g,norm_ratio=%.17g",
        stage, micro, orbital + 1, orbital < n_inactive ? "inactive" : "active", count,
        geometry.alpha, geometry.theta, geometry.step_norm,
        geometry.actual_angle, geometry.norm_ratio);
  }
}

}  // namespace xmvb::vb::diagnostics
