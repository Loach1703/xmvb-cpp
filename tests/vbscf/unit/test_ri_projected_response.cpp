#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "integrals/ao/ri_transform.hpp"
#include "integrals/response/ri_ao.hpp"

namespace {

using Matrix = Eigen::MatrixXd;
using xmvb::integrals::RiActiveTwoElectronResponseCache;
using xmvb::integrals::RiAoResponseWorkspace;
using xmvb::integrals::RiInactiveFactorView;

void require_close(const Matrix& actual, const Matrix& expected,
                   const std::string& name) {
  if (actual.rows() != expected.rows() || actual.cols() != expected.cols()) {
    throw std::runtime_error(name + " has the wrong shape");
  }
  if (actual.size() == 0) return;
  const double scale = std::max(1.0, expected.cwiseAbs().maxCoeff());
  const double error = (actual - expected).cwiseAbs().maxCoeff();
  if (error > 3.0e-12 * scale) {
    throw std::runtime_error(name + " differs from explicit AO oracle by " +
                             std::to_string(error));
  }
}

Matrix explicit_factor(const Matrix& packed, int auxiliary, int n) {
  Matrix factor = Matrix::Zero(n, n);
  for (int mu = 0; mu < n; ++mu) {
    for (int nu = 0; nu <= mu; ++nu) {
      const double value = packed(auxiliary, mu * (mu + 1) / 2 + nu);
      factor(mu, nu) = value;
      factor(nu, mu) = value;
    }
  }
  return factor;
}

Matrix explicit_ri_action(const Matrix& packed, const Matrix& density) {
  const int n = density.rows();
  Matrix result = Matrix::Zero(n, n);
  for (int auxiliary = 0; auxiliary < packed.rows(); ++auxiliary) {
    const Matrix factor = explicit_factor(packed, auxiliary, n);
    const double charge = (factor.array() * density.array()).sum();
    result.noalias() += 2.0 * charge * factor;
    result.noalias() -= factor * density * factor;
  }
  return result;
}

Matrix values(int rows, int cols, double phase) {
  Matrix result(rows, cols);
  for (Eigen::Index i = 0; i < result.size(); ++i) {
    result.data()[i] = 0.27 * std::sin(0.31 * (i + 1) + phase) +
                       0.13 * std::cos(0.17 * (i + 3) - phase);
  }
  return result;
}

Matrix pair_factors(const Matrix& packed, const Matrix& c) {
  const int na = c.cols();
  Matrix result(packed.rows(), na * (na + 1) / 2);
  for (int auxiliary = 0; auxiliary < packed.rows(); ++auxiliary) {
    const Matrix l = explicit_factor(packed, auxiliary, c.rows());
    for (int a = 0; a < na; ++a)
      for (int b = 0; b <= a; ++b)
        result(auxiliary, a * (a + 1) / 2 + b) =
            c.col(a).dot(l * c.col(b));
  }
  return result;
}

double explicit_pair_energy(const Matrix& packed, const Matrix& c,
                            const std::vector<double>& adjoint) {
  const Matrix factors = pair_factors(packed, c);
  double energy = 0.0;
  for (int p = 0; p < factors.cols(); ++p)
    for (int q = 0; q <= p; ++q)
      energy += adjoint[p * (p + 1) / 2 + q] *
          factors.col(p).dot(factors.col(q));
  return energy;
}

Matrix finite_difference_pair_gradient(const Matrix& packed, const Matrix& c,
                                       const std::vector<double>& adjoint) {
  constexpr double step = 1.0e-5;
  Matrix gradient(c.rows(), c.cols());
  for (int a = 0; a < c.cols(); ++a)
    for (int mu = 0; mu < c.rows(); ++mu) {
      Matrix plus = c, minus = c;
      plus(mu, a) += step;
      minus(mu, a) -= step;
      gradient(mu, a) =
          (explicit_pair_energy(packed, plus, adjoint) -
           explicit_pair_energy(packed, minus, adjoint)) / (2.0 * step);
    }
  return gradient;
}

void check_case(int full_count, int moving_count, bool cached,
                std::size_t budget) {
  constexpr int n = 6;
  constexpr int na = 2;
  constexpr int n_aux = 9;
  Matrix packed = values(n_aux, n * (n + 1) / 2, 0.12);
  Matrix c = values(n, na, 0.43);
  Matrix a0 = values(n, full_count, 0.71);
  Matrix x = values(n, moving_count, 1.27);
  Matrix w = values(n, na, -0.29);
  Matrix a(n, moving_count);
  const std::array<int, 2> map{2, 0};
  for (int moving = 0; moving < moving_count; ++moving) {
    a.col(moving) = a0.col(map[moving]);
  }

  RiActiveTwoElectronResponseCache cache;
  cache.n_basis_functions = n;
  cache.n_active_orbitals = na;
  cache.n_auxiliary_functions = n_aux;
  cache.metric_whitened_ao_pair_factors = &packed;
  cache.accepted_active_coefficients = &c;
  if (cached) {
    cache.accepted_one_index_factors = std::make_shared<Matrix>(
        xmvb::integrals::transform_ri_one_index_tiled(packed, c));
  }
  Matrix full_factors;
  RiInactiveFactorView inactive;
  if (cached && full_count != 0) {
    full_factors = xmvb::integrals::transform_ri_one_index_tiled(packed, a0);
    inactive.factors = &full_factors;
  }
  inactive.orbitals = std::span<const int>(map.data(), moving_count);

  const Matrix accepted_pair_factors = pair_factors(packed, c);
  cache.accepted_active_pair_factors = &accepted_pair_factors;
  for (int first = 0; first < na; ++first)
    for (int second = 0; second <= first; ++second) {
      cache.active_pair_first_indices.push_back(first);
      cache.active_pair_second_indices.push_back(second);
    }
  cache.auxiliary_tile_size = 3;

  const Matrix d = x * a.transpose() + a * x.transpose();
  const Matrix k = w * c.transpose() + c * w.transpose();
  const Matrix expected_inactive =
      2.0 * explicit_ri_action(packed, 2.0 * d + 0.5 * k) * a0;
  const Matrix expected_active = explicit_ri_action(packed, d) * c;

  RiAoResponseWorkspace workspace;
  const auto projected = xmvb::integrals::apply_ri_projected_ao_response(
      cache, a, x, w, a0, workspace, budget, inactive);
  require_close(projected.inactive, expected_inactive,
                "projected inactive RI response");
  require_close(projected.active, expected_active,
                "projected active RI response");
  // Reuse the same workspace with a different source and verify stale charges
  // or worker accumulators cannot leak into the next action.
  const Matrix zero_x = Matrix::Zero(n, moving_count);
  const Matrix zero_w = Matrix::Zero(n, na);
  const auto zero = xmvb::integrals::apply_ri_projected_ao_response(
      cache, a, zero_x, zero_w, a0, workspace, budget, inactive);
  require_close(zero.inactive, Matrix::Zero(n, full_count),
                "reused projected inactive RI response");
  require_close(zero.active, Matrix::Zero(n, na),
                "reused projected active RI response");

  if (full_count == 0) return;
  const Matrix active_direction = values(n, na, 0.93);
  const Matrix raw_weight = values(na, na, -0.37);
  const Matrix weight = 0.5 * (raw_weight + raw_weight.transpose()).eval();
  const Matrix fused_weighted = active_direction * weight;
  const std::vector<double> adjoint{0.14, -0.21, 0.09, 0.27, -0.12, 0.19};
  const xmvb::integrals::RiActiveHvpInput fixed{
      active_direction, weight, adjoint};
  const auto fused = xmvb::integrals::apply_ri_projected_ao_response(
      cache, a, x, fused_weighted, a0, workspace, budget, inactive, &fixed);
  const Matrix fused_k = fused_weighted * c.transpose() +
      c * fused_weighted.transpose();
  require_close(fused.inactive,
                2.0 * explicit_ri_action(packed, 2.0 * d + 0.5 * fused_k) * a0,
                "fused inactive RI response");
  require_close(fused.active, expected_active, "fused active RI response");
  require_close(fused.factor_direction,
                xmvb::integrals::compute_ri_active_pair_factor_directional_derivative(
                    cache, active_direction), "fused active-pair factor direction");
  require_close(fused.fixed_two_electron,
                xmvb::integrals::apply_ri_packed_active_two_electron_adjoint_hessian_vector(
                    adjoint, cache, active_direction), "fused fixed 2e image");
  if (cached) {
    RiInactiveFactorView generated_inactive = inactive;
    generated_inactive.factors = nullptr;
    const auto mixed = xmvb::integrals::apply_ri_projected_ao_response(
        cache, a, x, fused_weighted, a0, workspace, budget,
        generated_inactive, &fixed);
    require_close(mixed.inactive, fused.inactive,
                  "resident-active/generated-inactive fused response");
    require_close(mixed.active, fused.active,
                  "resident-active/generated-inactive active response");
    require_close(mixed.factor_direction, fused.factor_direction,
                  "resident-active/generated-inactive factor direction");
    require_close(mixed.fixed_two_electron, fused.fixed_two_electron,
                  "resident-active/generated-inactive fixed image");
  }
  constexpr double step = 1.0e-4;
  const Matrix finite_difference =
      (finite_difference_pair_gradient(packed, c + step * active_direction, adjoint) -
       finite_difference_pair_gradient(packed, c - step * active_direction, adjoint)) /
      (2.0 * step);
  if ((fused.fixed_two_electron - finite_difference).norm() >
      3.0e-4 * std::max(1.0, finite_difference.norm()))
    throw std::runtime_error("fused fixed 2e image disagrees with independent energy FD");

  const Matrix zero_direction = Matrix::Zero(n, na);
  const xmvb::integrals::RiActiveHvpInput zero_fixed{
      zero_direction, weight, adjoint};
  const auto zero_fused = xmvb::integrals::apply_ri_projected_ao_response(
      cache, a, x, zero_direction, a0, workspace, budget, inactive, &zero_fixed);
  require_close(zero_fused.factor_direction, Matrix::Zero(n_aux, na * (na + 1) / 2),
                "zero fused factor direction");
  require_close(zero_fused.fixed_two_electron, Matrix::Zero(n, na),
                "zero fused fixed 2e image");
}

}  // namespace

int main() {
#ifdef _OPENMP
  const int saved_threads = omp_get_max_threads();
#endif
  try {
    for (int threads : {1, 2, 4, 16}) {
#ifdef _OPENMP
      omp_set_num_threads(threads);
#else
      (void)threads;
#endif
      for (std::size_t budget : {std::size_t{8} * 1024,
                                 std::size_t{64} * 1024 * 1024}) {
        for (bool cached : {false, true}) {
          check_case(3, 2, cached, budget);
          check_case(3, 0, cached, budget);
          check_case(5, 1, cached, budget);
          check_case(5, 0, cached, budget);
          check_case(0, 0, cached, budget);
        }
      }
    }
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
#ifdef _OPENMP
  omp_set_num_threads(saved_threads);
#endif
  return 0;
}
