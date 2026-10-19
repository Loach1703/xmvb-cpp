#include "dft/integration/potential.hpp"

#include "core/workspace.hpp"
#include "dft/grid/atomic.hpp"
#include "dft/integration/density.hpp"
#include "dft/xc/functional.hpp"
#include "integrals/basis/values.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>

namespace xmvb::dft {
namespace {

void compensated_add(double term, double* sum, double* correction) {
  const double adjusted = term - *correction;
  const double next = *sum + adjusted;
  *correction = (next - *sum) - adjusted;
  *sum = next;
}

}  // namespace

LdaIntegrator::LdaIntegrator(const integrals::LibcintInput& input,
                             const Grid& grid,
                             const LdaFunctional& functional)
    : grid_(grid), functional_(functional),
      basis_(std::make_unique<integrals::BasisEvaluator>(input)) {
  if (grid_.points.rows() == 0 || grid_.points.cols() != 3 ||
      grid_.weights.size() != grid_.points.rows() ||
      !grid_.points.allFinite() || !grid_.weights.allFinite() ||
      (grid_.weights.array() < 0.0).any())
    throw std::invalid_argument("invalid DFT grid geometry or weights");
}

LdaIntegrator::~LdaIntegrator() = default;

XcResult LdaIntegrator::evaluate(
    const Eigen::Ref<const Eigen::MatrixXd>& density) const {
  const Eigen::Index n = basis_->n_bf();
  if (density.rows() != n || density.cols() != n || !density.allFinite() ||
      !density.isApprox(density.transpose(), 1.0e-12)) {
    throw std::invalid_argument("restricted AO density must be finite and symmetric");
  }

  const std::size_t n_basis = static_cast<std::size_t>(n);
  const std::size_t max_bytes = std::numeric_limits<std::size_t>::max();
  if (n_basis && n_basis > max_bytes / sizeof(double) / n_basis) {
    throw std::length_error("DFT XC matrix size overflow");
  }
  const std::size_t output_bytes = n_basis * n_basis * sizeof(double);
  if (n_basis > (max_bytes / sizeof(double) - 6) / 2) {
    throw std::length_error("DFT XC tile size overflow");
  }
  // Per point: AO/projected AO (2*n), rho/exc/vrho (3), and Libxc input
  // plus two component arrays (3). Coordinates are a borrowed block view.
  const std::size_t point_bytes = (2 * n_basis + 6) * sizeof(double);
  const std::size_t available = workspace::available();
  if (output_bytes > available ||
      point_bytes > available - output_bytes) {
    throw std::runtime_error("DFT XC resource_limit: insufficient shared workspace");
  }
  const std::size_t grid_size =
      static_cast<std::size_t>(grid_.points.rows());
  const std::size_t rows_per_tile =
      std::min(grid_size, (available - output_bytes) / point_bytes);
  workspace::Reservation reservation;
  if (!reservation.try_resize(output_bytes + rows_per_tile * point_bytes)) {
    throw std::runtime_error("DFT XC resource_limit: workspace reservation failed");
  }

  XcResult result;
  result.potential = Eigen::MatrixXd::Zero(n, n);
  integrals::BasisValues basis_values;
  Eigen::MatrixXd projected;
  Eigen::VectorXd rho;
  Eigen::VectorXd exc;
  Eigen::VectorXd vrho;
  double energy_correction = 0.0;
  double electron_correction = 0.0;
  for (Eigen::Index begin = 0; begin < grid_.points.rows();) {
    const Eigen::Index rows = std::min<Eigen::Index>(
        grid_.points.rows() - begin,
        static_cast<Eigen::Index>(rows_per_tile));
    basis_->evaluate(grid_.points.middleRows(begin, rows), 0, &basis_values);
    const auto& ao = basis_values.value;
    evaluate_density_tile(ao, density, &projected, &rho);
    functional_.evaluate(rho, &exc, &vrho);

    for (Eigen::Index point = 0; point < rows; ++point) {
      const double weight = grid_.weights[begin + point];
      compensated_add(weight * rho[point], &result.electron_count,
                      &electron_correction);
      compensated_add(weight * rho[point] * exc[point], &result.energy,
                      &energy_correction);
      const double weighted_vrho = weight * vrho[point];
      if (!std::isfinite(weighted_vrho)) {
        throw std::runtime_error("nonfinite weighted LDA potential");
      }
      projected.row(point) = weighted_vrho * ao.row(point);
    }
    result.potential.noalias() += ao.transpose() * projected;
    begin += rows;
  }
  if (!std::isfinite(result.energy) ||
      !std::isfinite(result.electron_count) ||
      !result.potential.allFinite()) {
    throw std::runtime_error("nonfinite LDA grid integral");
  }
  // Remove only roundoff asymmetry in the symmetric AO contraction.
  for (Eigen::Index i = 0; i < n; ++i) {
    for (Eigen::Index j = 0; j < i; ++j) {
      const double value = 0.5 * (result.potential(i, j) +
                                  result.potential(j, i));
      result.potential(i, j) = value;
      result.potential(j, i) = value;
    }
  }
  return result;
}

}  // namespace xmvb::dft
