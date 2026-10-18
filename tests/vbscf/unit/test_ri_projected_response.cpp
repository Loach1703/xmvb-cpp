#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>

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
}

}  // namespace

int main() {
#ifdef _OPENMP
  const int saved_threads = omp_get_max_threads();
#endif
  try {
    for (int threads : {1, 2, 4}) {
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
