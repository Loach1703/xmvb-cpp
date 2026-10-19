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

void symmetrize_sum(Eigen::MatrixXd* matrix) {
  for (Eigen::Index i = 0; i < matrix->rows(); ++i) {
    (*matrix)(i, i) *= 2.0;
    for (Eigen::Index j = 0; j < i; ++j) {
      const double value = (*matrix)(i, j) + (*matrix)(j, i);
      (*matrix)(i, j) = value;
      (*matrix)(j, i) = value;
    }
  }
}

}  // namespace

XcIntegrator::XcIntegrator(const integrals::LibcintInput& input,
                           const Grid& grid,
                           const Functional& functional)
    : grid_(grid), functional_(functional),
      basis_(std::make_unique<integrals::BasisEvaluator>(input)) {
  if (grid_.points.rows() == 0 || grid_.points.cols() != 3 ||
      grid_.weights.size() != grid_.points.rows() ||
      !grid_.points.allFinite() || !grid_.weights.allFinite() ||
      (grid_.weights.array() < 0.0).any()) {
    throw std::invalid_argument("invalid DFT grid geometry or weights");
  }
}

XcIntegrator::~XcIntegrator() = default;

XcResult XcIntegrator::evaluate(
    std::span<const Eigen::Ref<const Eigen::MatrixXd>> density) const {
  const int spins = functional_.n_spin();
  const bool gga = functional_.is_gga();
  const int sigma_columns = gga ? (spins == 1 ? 1 : 3) : 0;
  const Eigen::Index n = basis_->n_bf();
  if (density.size() != static_cast<std::size_t>(spins)) {
    throw std::invalid_argument("AO density spin count does not match XC functional");
  }
  for (const auto& matrix : density) {
    if (matrix.rows() != n || matrix.cols() != n || !matrix.allFinite() ||
        !matrix.isApprox(matrix.transpose(), 1.0e-12)) {
      throw std::invalid_argument("AO density must be finite and symmetric");
    }
  }

  const std::size_t nb = static_cast<std::size_t>(n);
  const std::size_t max_doubles =
      std::numeric_limits<std::size_t>::max() / sizeof(double);
  if (nb && nb > max_doubles / nb / static_cast<std::size_t>(spins)) {
    throw std::length_error("DFT XC matrix size overflow");
  }
  const std::size_t output_doubles =
      static_cast<std::size_t>(spins) * nb * nb;
  // Simultaneously live tile arrays:
  // LDA: A, projected/AP, rho, exc, vrho, Libxc's packed rho/exc/vrho.
  // GGA: A+3 grad(A), projected/AP, rho+3 grad(rho), sigma,
  //      exc+vrho+vsigma, and Libxc's packed rho/sigma/exc/vrho/vsigma.
  const std::size_t basis_factor = gga ? 5 : 2;
  const std::size_t scalar_factor = gga
      ? static_cast<std::size_t>(2 + 7 * spins + 4 * sigma_columns)
      : static_cast<std::size_t>(2 + 4 * spins);
  if (nb > (max_doubles - scalar_factor) / basis_factor) {
    throw std::length_error("DFT XC tile size overflow");
  }
  const std::size_t point_doubles = basis_factor * nb + scalar_factor;
  const std::size_t available = workspace::available();
  if (output_doubles > available / sizeof(double) ||
      point_doubles > (available / sizeof(double) - output_doubles)) {
    throw std::runtime_error("DFT XC resource_limit: insufficient shared workspace");
  }
  const std::size_t grid_size =
      static_cast<std::size_t>(grid_.points.rows());
  const std::size_t rows_per_tile = std::min(
      grid_size, (available / sizeof(double) - output_doubles) / point_doubles);
  workspace::Reservation reservation;
  if (!reservation.try_resize(
          (output_doubles + rows_per_tile * point_doubles) * sizeof(double))) {
    throw std::runtime_error("DFT XC resource_limit: workspace reservation failed");
  }

  XcResult result;
  for (int spin = 0; spin < spins; ++spin) {
    result.potential[spin] = Eigen::MatrixXd::Zero(n, n);
  }
  integrals::BasisValues basis_values;
  Eigen::MatrixXd projected;
  Eigen::MatrixXd rho;
  std::array<Eigen::MatrixXd, 3> rho_gradient;
  Eigen::MatrixXd sigma;
  XcValues xc;
  double energy_correction = 0.0;
  std::array<double, 2> electron_correction{};

  for (Eigen::Index begin = 0; begin < grid_.points.rows();) {
    const Eigen::Index rows = std::min<Eigen::Index>(
        grid_.points.rows() - begin,
        static_cast<Eigen::Index>(rows_per_tile));
    basis_->evaluate(grid_.points.middleRows(begin, rows), gga ? 1 : 0,
                     &basis_values);
    const auto& ao = basis_values.value;
    rho.resize(rows, spins);
    if (gga) {
      for (auto& component : rho_gradient) component.resize(rows, spins);
    }
    for (int spin = 0; spin < spins; ++spin) {
      evaluate_density_tile(basis_values, density[spin], gga, spin,
                            &projected, &rho, &rho_gradient);
    }

    sigma.resize(rows, sigma_columns);
    if (gga) {
      for (Eigen::Index point = 0; point < rows; ++point) {
        double aa = 0.0;
        double ab = 0.0;
        double bb = 0.0;
        for (int axis = 0; axis < 3; ++axis) {
          const double alpha = rho_gradient[axis](point, 0);
          aa += alpha * alpha;
          if (spins == 2) {
            const double beta = rho_gradient[axis](point, 1);
            ab += alpha * beta;
            bb += beta * beta;
          }
        }
        sigma(point, 0) = aa;
        if (spins == 2) {
          sigma(point, 1) = ab;
          sigma(point, 2) = bb;
        }
      }
      if (!sigma.allFinite()) {
        throw std::runtime_error("nonfinite density-gradient invariant on DFT grid");
      }
    }
    functional_.evaluate(rho, sigma, &xc);

    for (Eigen::Index point = 0; point < rows; ++point) {
      const double weight = grid_.weights[begin + point];
      double total_rho = 0.0;
      for (int spin = 0; spin < spins; ++spin) {
        const double spin_rho = rho(point, spin);
        total_rho += spin_rho;
        compensated_add(weight * spin_rho, &result.electron_count[spin],
                        &electron_correction[spin]);
      }
      compensated_add(weight * total_rho * xc.exc[point], &result.energy,
                      &energy_correction);
    }

    // projected was used for A*P above. Reuse it for the weak-form
    // contribution and contract each spin with a single AO transpose GEMM.
    for (int spin = 0; spin < spins; ++spin) {
      for (Eigen::Index point = 0; point < rows; ++point) {
        const double weight = grid_.weights[begin + point];
        const double density_scale = 0.5 * weight * xc.vrho(point, spin);
        std::array<double, 3> gradient_scale{};
        if (gga) {
          for (int axis = 0; axis < 3; ++axis) {
            const double own = rho_gradient[axis](point, spin);
            double t = 0.0;
            if (spins == 1) {
              t = 2.0 * xc.vsigma(point, 0) * own;
            } else if (spin == 0) {
              t = 2.0 * xc.vsigma(point, 0) * own +
                  xc.vsigma(point, 1) * rho_gradient[axis](point, 1);
            } else {
              t = 2.0 * xc.vsigma(point, 2) * own +
                  xc.vsigma(point, 1) * rho_gradient[axis](point, 0);
            }
            gradient_scale[axis] = weight * t;
          }
        }
        for (Eigen::Index orbital = 0; orbital < n; ++orbital) {
          double value = density_scale * ao(point, orbital);
          if (gga) {
            for (int axis = 0; axis < 3; ++axis) {
              value += gradient_scale[axis] *
                       basis_values.gradient[axis](point, orbital);
            }
          }
          if (!std::isfinite(value)) {
            throw std::runtime_error("nonfinite weighted XC potential");
          }
          projected(point, orbital) = value;
        }
      }
      result.potential[spin].noalias() += ao.transpose() * projected;
    }
    begin += rows;
  }

  if (!std::isfinite(result.energy)) {
    throw std::runtime_error("nonfinite XC grid integral");
  }
  for (int spin = 0; spin < spins; ++spin) {
    symmetrize_sum(&result.potential[spin]);
    if (!std::isfinite(result.electron_count[spin]) ||
        !result.potential[spin].allFinite()) {
      throw std::runtime_error("nonfinite XC grid integral");
    }
  }
  return result;
}

std::array<Eigen::MatrixXd, 2> XcIntegrator::response(
    std::span<const Eigen::Ref<const Eigen::MatrixXd>> density,
    std::span<const Eigen::Ref<const Eigen::MatrixXd>> delta_density) const {
  const int spins = functional_.n_spin();
  const bool gga = functional_.is_gga();
  const int sigma_columns = gga ? (spins == 1 ? 1 : 3) : 0;
  const Eigen::Index n = basis_->n_bf();
  if (density.size() != static_cast<std::size_t>(spins) ||
      delta_density.size() != static_cast<std::size_t>(spins)) {
    throw std::invalid_argument("AO density spin count does not match XC functional");
  }
  const auto validate_matrix = [n](
      const Eigen::Ref<const Eigen::MatrixXd>& matrix) {
    if (matrix.rows() != n || matrix.cols() != n || !matrix.allFinite() ||
        !matrix.isApprox(matrix.transpose(), 1.0e-12)) {
      throw std::invalid_argument(
          "AO density and direction must be finite and symmetric");
    }
  };
  for (int spin = 0; spin < spins; ++spin) {
    validate_matrix(density[spin]);
    validate_matrix(delta_density[spin]);
  }

  const std::size_t nb = static_cast<std::size_t>(n);
  const std::size_t max_doubles =
      std::numeric_limits<std::size_t>::max() / sizeof(double);
  if (nb && nb > max_doubles / nb / static_cast<std::size_t>(spins)) {
    throw std::length_error("DFT XC matrix size overflow");
  }
  const std::size_t output_doubles =
      static_cast<std::size_t>(spins) * nb * nb;
  // Basis values, gradients and the one reusable projected AO matrix are
  // the only point-by-AO arrays. The scalar allowance covers base and
  // directional densities/invariants, XC outputs, and Libxc buffers.
  const std::size_t basis_factor = gga ? 5 : 2;
  const std::size_t scalar_factor = gga
      ? static_cast<std::size_t>(2 + 20 * spins + 14 * sigma_columns)
      : static_cast<std::size_t>(2 + 10 * spins);
  if (nb > (max_doubles - scalar_factor) / basis_factor) {
    throw std::length_error("DFT XC tile size overflow");
  }
  const std::size_t point_doubles = basis_factor * nb + scalar_factor;
  const std::size_t available = workspace::available();
  if (output_doubles > available / sizeof(double) ||
      point_doubles > (available / sizeof(double) - output_doubles)) {
    throw std::runtime_error("DFT XC response resource_limit: insufficient shared workspace");
  }
  const std::size_t grid_size =
      static_cast<std::size_t>(grid_.points.rows());
  const std::size_t rows_per_tile = std::min(
      grid_size, (available / sizeof(double) - output_doubles) / point_doubles);
  workspace::Reservation reservation;
  if (!reservation.try_resize(
          (output_doubles + rows_per_tile * point_doubles) * sizeof(double))) {
    throw std::runtime_error("DFT XC response resource_limit: workspace reservation failed");
  }

  std::array<Eigen::MatrixXd, 2> result;
  for (int spin = 0; spin < spins; ++spin) {
    result[spin] = Eigen::MatrixXd::Zero(n, n);
  }
  integrals::BasisValues basis_values;
  Eigen::MatrixXd projected;
  Eigen::MatrixXd rho, delta_rho, sigma, delta_sigma;
  std::array<Eigen::MatrixXd, 3> rho_gradient, delta_rho_gradient;
  XcValues xc;
  Eigen::MatrixXd delta_vrho, delta_vsigma;

  for (Eigen::Index begin = 0; begin < grid_.points.rows();) {
    const Eigen::Index rows = std::min<Eigen::Index>(
        grid_.points.rows() - begin,
        static_cast<Eigen::Index>(rows_per_tile));
    basis_->evaluate(grid_.points.middleRows(begin, rows), gga ? 1 : 0,
                     &basis_values);
    const auto& ao = basis_values.value;
    rho.resize(rows, spins);
    delta_rho.resize(rows, spins);
    if (gga) {
      for (int axis = 0; axis < 3; ++axis) {
        rho_gradient[axis].resize(rows, spins);
        delta_rho_gradient[axis].resize(rows, spins);
      }
    }
    for (int spin = 0; spin < spins; ++spin) {
      evaluate_density_tile(basis_values, density[spin], gga, spin,
                            &projected, &rho, &rho_gradient);
      evaluate_density_tile(basis_values, delta_density[spin], gga, spin,
                            &projected, &delta_rho, &delta_rho_gradient,
                            DensityTileMode::Directional);
    }

    sigma.resize(rows, sigma_columns);
    delta_sigma.resize(rows, sigma_columns);
    if (gga) {
      for (Eigen::Index point = 0; point < rows; ++point) {
        double aa = 0.0, ab = 0.0, bb = 0.0;
        double delta_aa = 0.0, delta_ab = 0.0, delta_bb = 0.0;
        for (int axis = 0; axis < 3; ++axis) {
          const double alpha = rho_gradient[axis](point, 0);
          const double delta_alpha = delta_rho_gradient[axis](point, 0);
          aa += alpha * alpha;
          delta_aa += 2.0 * alpha * delta_alpha;
          if (spins == 2) {
            const double beta = rho_gradient[axis](point, 1);
            const double delta_beta = delta_rho_gradient[axis](point, 1);
            ab += alpha * beta;
            delta_ab += delta_alpha * beta + alpha * delta_beta;
            bb += beta * beta;
            delta_bb += 2.0 * beta * delta_beta;
          }
        }
        sigma(point, 0) = aa;
        delta_sigma(point, 0) = delta_aa;
        if (spins == 2) {
          sigma(point, 1) = ab;
          sigma(point, 2) = bb;
          delta_sigma(point, 1) = delta_ab;
          delta_sigma(point, 2) = delta_bb;
        }
      }
      if (!sigma.allFinite() || !delta_sigma.allFinite()) {
        throw std::runtime_error("nonfinite density-gradient invariant on DFT grid");
      }
    }
    if (gga) functional_.evaluate(rho, sigma, &xc);
    functional_.response(rho, sigma, delta_rho, delta_sigma,
                         &delta_vrho, &delta_vsigma);

    // Reuse projected for the differentiated weak form. The final AO
    // transpose contraction and symmetrization match evaluate().
    for (int spin = 0; spin < spins; ++spin) {
      for (Eigen::Index point = 0; point < rows; ++point) {
        const double weight = grid_.weights[begin + point];
        const double density_scale = 0.5 * weight * delta_vrho(point, spin);
        std::array<double, 3> gradient_scale{};
        if (gga) {
          const int own_channel = spins == 1 ? 0 : (spin == 0 ? 0 : 2);
          for (int axis = 0; axis < 3; ++axis) {
            const double own = rho_gradient[axis](point, spin);
            const double delta_own = delta_rho_gradient[axis](point, spin);
            double t = 2.0 * (delta_vsigma(point, own_channel) * own +
                              xc.vsigma(point, own_channel) * delta_own);
            if (spins == 2) {
              const int other = 1 - spin;
              t += delta_vsigma(point, 1) *
                       rho_gradient[axis](point, other) +
                   xc.vsigma(point, 1) *
                       delta_rho_gradient[axis](point, other);
            }
            gradient_scale[axis] = weight * t;
          }
        }
        for (Eigen::Index orbital = 0; orbital < n; ++orbital) {
          double value = density_scale * ao(point, orbital);
          if (gga) {
            for (int axis = 0; axis < 3; ++axis) {
              value += gradient_scale[axis] *
                       basis_values.gradient[axis](point, orbital);
            }
          }
          if (!std::isfinite(value)) {
            throw std::runtime_error("nonfinite weighted XC response");
          }
          projected(point, orbital) = value;
        }
      }
      result[spin].noalias() += ao.transpose() * projected;
    }
    begin += rows;
  }

  for (int spin = 0; spin < spins; ++spin) {
    symmetrize_sum(&result[spin]);
    if (!result[spin].allFinite()) {
      throw std::runtime_error("nonfinite XC grid response");
    }
  }
  return result;
}

}  // namespace xmvb::dft
