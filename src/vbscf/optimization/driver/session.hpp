#pragma once

#include <Eigen/Core>

#include "vbscf/optimization/driver/options.hpp"
#include "vbscf/orbitals/charts/chart.hpp"
#include "vbscf/optimization/preconditioners/type.hpp"

namespace xmvb::vb {

class SparseParameterLayout;
class VbScfObjective;
struct VbScfOptimizerResult;
struct TnhvpIterationRecord;

namespace optimizer_detail {

OrbitalChart build_orbital_chart(
    const VbScfObjective& objective,
    const SparseParameterLayout& parameter_view,
    OrbitalPreconditioner preconditioner = OrbitalPreconditioner::OneElectron,
    bool need_ah_mapping = false);

/**
 * @brief Certifies the shared stopping target in the accepted orbital chart.
 *
 * Complete OEO also requires the gauge-invariant quotient dual norm. Other
 * manifolds retain their established projected-gradient infinity criterion.
 */
bool orbital_gradient_converged(
    const VbScfObjective& objective,
    const SparseParameterLayout& parameter_view,
    const OrbitalChart& chart,
    const Eigen::VectorXd& reduced_gradient,
    double tolerance);

void sync_result_from_objective(
    const VbScfObjective& objective,
    VbScfOptimizerResult* result);

void record_accepted_iteration_snapshot(
    VbScfObjective* objective,
    int accepted_iteration_index,
    const VbScfOptimizerOptions& options,
    VbScfOptimizerResult* result,
    const TnhvpIterationRecord* tnhvp = nullptr,
    const Eigen::VectorXd* reduced_gradient = nullptr);

}  // namespace optimizer_detail

}  // namespace xmvb::vb
