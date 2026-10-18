#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include "core/pair_index.hpp"
#include "integrals/active/exact.hpp"
#include "integrals/active/ri.hpp"
#include "integrals/active/two_electron.hpp"
#include "integrals/active/transform.hpp"
#include "integrals/ao/graph.hpp"
#include "integrals/response/adjoint.hpp"
#include "integrals/response/ri.hpp"
#include "vbscf/structures/orthogonal_ci/integrals.hpp"

namespace {

using xmvb::integrals::ActiveSpaceTwoElectronResult;
using xmvb::integrals::AoIntegralInput;
using xmvb::integrals::RiAoFactorization;

void require_close(const Eigen::MatrixXd& actual,
    const Eigen::MatrixXd& expected, double tolerance, const char* label) {
  if (actual.rows() != expected.rows() || actual.cols() != expected.cols() ||
      (actual - expected).norm() >
          tolerance * std::max(1.0, expected.norm())) {
    throw std::runtime_error(label);
  }
}

std::vector<double> pack_pair_adjoint(const Eigen::MatrixXd& gradient) {
  const int pairs = gradient.rows();
  std::vector<double> packed(pairs * (pairs + 1) / 2);
  for (int column = 0; column < pairs; ++column)
    for (int row = 0; row <= column; ++row)
      packed[xmvb::core::PairIndex::packed_pair_of_pairs_index(row, column)] =
          row == column ? gradient(row, column)
                        : gradient(row, column) + gradient(column, row);
  return packed;
}

ActiveSpaceTwoElectronResult exact_at(
    const AoIntegralInput& ao, const Eigen::MatrixXd& coefficients) {
  const int n_basis = coefficients.rows();
  const int n_active = coefficients.cols();
  std::vector<int> row_offsets(n_basis + 1, 0);
  std::vector<int> orbital_indices;
  std::vector<double> values;
  for (int row = 0; row < n_basis; ++row) {
    for (int column = 0; column < n_active; ++column) {
      const double value = coefficients(row, column);
      if (value != 0.0) {
        orbital_indices.push_back(column);
        values.push_back(value);
      }
    }
    row_offsets[row + 1] = static_cast<int>(values.size());
  }
  return xmvb::integrals::transform_exact_active_integrals(
      ao, row_offsets, orbital_indices, values, n_active);
}

void check_frame(bool ri_backend) {
  constexpr int n_basis = 4;
  constexpr int n_active = 3;
  constexpr int n_aux = 5;
  constexpr int ao_pairs = n_basis * (n_basis + 1) / 2;
  constexpr int active_pairs = n_active * (n_active + 1) / 2;
  RiAoFactorization ri;
  ri.n_basis_functions = n_basis;
  ri.n_auxiliary_functions = n_aux;
  ri.n_packed_ao_pairs = ao_pairs;
  ri.metric_whitened_ao_pair_factors.resize(n_aux, ao_pairs);
  for (int column = 0; column < ao_pairs; ++column)
    for (int row = 0; row < n_aux; ++row)
      ri.metric_whitened_ao_pair_factors(row, column) =
          std::sin(0.19 * (row + 1) * (column + 2));
  const Eigen::MatrixXd ao_kernel =
      ri.metric_whitened_ao_pair_factors.transpose() *
      ri.metric_whitened_ao_pair_factors;
  std::vector<int> left, right;
  std::vector<double> values;
  for (int p = 0; p < ao_pairs; ++p)
    for (int q = 0; q <= p; ++q) {
      left.push_back(p);
      right.push_back(q);
      values.push_back(ao_kernel(p, q));
    }
  AoIntegralInput ao;
  ao.n_basis_functions = n_basis;
  ao.pair_graph = xmvb::integrals::build_ao_pair_graph(
      left, right, values, n_basis);

  Eigen::MatrixXd c(n_basis, n_active), dc(n_basis, n_active);
  c << 1.1, 0.2, 0.0,
       0.3, 1.0, -0.1,
       0.0, 0.2, 0.9,
       -0.2, 0.0, 0.4;
  for (int column = 0; column < n_active; ++column)
    for (int row = 0; row < n_basis; ++row)
      dc(row, column) = 0.04 * std::cos(0.31 * (row + 2) * (column + 1));
  Eigen::MatrixXd ao_overlap = Eigen::MatrixXd::Identity(n_basis, n_basis);
  for (int column = 0; column < n_basis; ++column)
    for (int row = 0; row < n_basis; ++row)
      ao_overlap(row, column) += 0.03 * std::cos(0.2 * (row + 1) * (column + 1));
  const Eigen::MatrixXd s = c.transpose() * ao_overlap * c;
  const Eigen::MatrixXd ds =
      dc.transpose() * ao_overlap * c + c.transpose() * ao_overlap * dc;
  const Eigen::MatrixXd r = s.llt().matrixU();
  const Eigen::MatrixXd q = r.triangularView<Eigen::Upper>().solve(
      Eigen::MatrixXd::Identity(n_active, n_active));
  const Eigen::MatrixXd cbar = c * q;
  Eigen::MatrixXd h(n_active, n_active), dh(n_active, n_active);
  for (int column = 0; column < n_active; ++column)
    for (int row = 0; row < n_active; ++row) {
      h(row, column) = 0.2 * std::sin(0.3 * (row + column + 2));
      dh(row, column) = 0.01 * std::cos(0.4 * (row + column + 2));
    }
  const auto forward = [&](const Eigen::MatrixXd& coefficients) {
    return ri_backend
        ? xmvb::integrals::transform_ri_active_integrals(ri, coefficients)
        : exact_at(ao, coefficients);
  };
  const auto raw = forward(c);
  auto framed = forward(cbar);
  framed.orbital_transform = r;
  const std::vector<double> overlap(s.data(), s.data() + s.size());
  const std::vector<double> overlap_direction(ds.data(), ds.data() + ds.size());
  const auto raw_orth = xmvb::vb::orthogonalize_active_integrals(
      overlap, h, raw, n_active);
  const auto framed_orth = xmvb::vb::orthogonalize_active_integrals(
      overlap, h, framed, n_active);
  require_close(framed_orth.one_electron, raw_orth.one_electron, 1e-12,
      "framed one-electron integrals differ");
  require_close(framed_orth.pair_kernel, raw_orth.pair_kernel, 2e-11,
      "framed two-electron integrals differ");
  if (ri_backend)
    require_close(framed_orth.ri_pair_factors, raw_orth.ri_pair_factors,
        2e-11, "framed RI factors differ");

  Eigen::MatrixXd one_adjoint(n_active, n_active),
      pair_adjoint(active_pairs, active_pairs), generator(n_active, n_active);
  for (int column = 0; column < n_active; ++column)
    for (int row = 0; row < n_active; ++row) {
      one_adjoint(row, column) = std::sin(0.13 * (row + 2) * (column + 1));
      generator(row, column) = std::cos(0.17 * (row + 1) * (column + 3));
    }
  for (int column = 0; column < active_pairs; ++column)
    for (int row = 0; row < active_pairs; ++row)
      pair_adjoint(row, column) =
          std::cos(0.11 * (row + 1) * (column + 2));
  const auto raw_adjoint = xmvb::vb::backpropagate_orthogonal_active_integral_adjoint(
      raw_orth, one_adjoint, pair_adjoint, generator);
  const auto framed_adjoint = xmvb::vb::backpropagate_orthogonal_active_integral_adjoint(
      framed_orth, one_adjoint, pair_adjoint, generator);
  require_close(framed_adjoint.overlap, raw_adjoint.overlap, 1e-12,
      "framed overlap adjoint differs");
  require_close(framed_adjoint.one_electron, raw_adjoint.one_electron, 1e-12,
      "framed one-electron adjoint differs");
  require_close(framed_adjoint.pair_kernel,
      0.5 * (pair_adjoint + pair_adjoint.transpose()), 1e-12,
      "framed pair adjoint was transformed twice");

  Eigen::MatrixXd raw_coefficient_gradient, framed_coefficient_gradient;
  if (ri_backend) {
    const auto raw_cache = xmvb::integrals::build_ri_active_two_electron_response_cache(
        ri, raw, n_active);
    const auto framed_cache = xmvb::integrals::build_ri_active_two_electron_response_cache(
        ri, framed, n_active);
    raw_coefficient_gradient = xmvb::integrals::backpropagate_ri_packed_active_two_electron_gradient(
        pack_pair_adjoint(raw_adjoint.pair_kernel), raw_cache);
    framed_coefficient_gradient = xmvb::integrals::backpropagate_ri_packed_active_two_electron_gradient(
        pack_pair_adjoint(framed_adjoint.pair_kernel), framed_cache);
  } else {
    raw_coefficient_gradient = xmvb::integrals::backpropagate_exact_packed_active_two_electron_gradient(
        pack_pair_adjoint(raw_adjoint.pair_kernel), ao, raw);
    framed_coefficient_gradient = xmvb::integrals::backpropagate_exact_packed_active_two_electron_gradient(
        pack_pair_adjoint(framed_adjoint.pair_kernel), ao, framed);
  }
  require_close(framed_coefficient_gradient * q.transpose(),
      raw_coefficient_gradient, 2e-10,
      "framed physical orbital pullback differs");

  constexpr double step = 1e-6;
  const auto raw_plus = forward(c + step * dc);
  const auto raw_minus = forward(c - step * dc);
  const auto framed_plus = forward(cbar + step * dc * q);
  const auto framed_minus = forward(cbar - step * dc * q);
  std::vector<double> raw_delta, framed_delta;
  Eigen::MatrixXd raw_factor_delta, framed_factor_delta;
  if (ri_backend) {
    raw_factor_delta = (raw_plus.ri_active_pair_factors -
        raw_minus.ri_active_pair_factors) / (2 * step);
    framed_factor_delta = (framed_plus.ri_active_pair_factors -
        framed_minus.ri_active_pair_factors) / (2 * step);
  } else {
    raw_delta.resize(raw_plus.packed_active_two_electron_integrals.size());
    framed_delta.resize(framed_plus.packed_active_two_electron_integrals.size());
    for (std::size_t index = 0; index < raw_delta.size(); ++index) {
      raw_delta[index] = (raw_plus.packed_active_two_electron_integrals[index] -
          raw_minus.packed_active_two_electron_integrals[index]) / (2 * step);
      framed_delta[index] = (framed_plus.packed_active_two_electron_integrals[index] -
          framed_minus.packed_active_two_electron_integrals[index]) / (2 * step);
    }
  }
  const auto raw_jet = xmvb::vb::orthogonalize_active_integral_direction(
      raw_orth, overlap_direction, dh, raw_delta, n_active,
      ri_backend ? &raw_factor_delta : nullptr);
  const auto framed_jet = xmvb::vb::orthogonalize_active_integral_direction(
      framed_orth, overlap_direction, dh, framed_delta, n_active,
      ri_backend ? &framed_factor_delta : nullptr);
  require_close(framed_jet.orbital_transform, raw_jet.orbital_transform, 1e-12,
      "framed Cholesky direction differs");
  require_close(framed_jet.one_electron, raw_jet.one_electron, 1e-12,
      "framed one-electron direction differs");
  require_close(framed_jet.pair_kernel, raw_jet.pair_kernel, 2e-8,
      "framed two-electron direction differs");

  const Eigen::MatrixXd zero_one = Eigen::MatrixXd::Zero(n_active, n_active);
  const Eigen::MatrixXd zero_pair = Eigen::MatrixXd::Zero(active_pairs, active_pairs);
  const auto raw_adjoint_jet = xmvb::vb::backpropagate_orthogonal_active_integral_adjoint_direction(
      raw_orth, raw_jet, one_adjoint, pair_adjoint, generator,
      zero_one, zero_pair, zero_one);
  const auto framed_adjoint_jet = xmvb::vb::backpropagate_orthogonal_active_integral_adjoint_direction(
      framed_orth, framed_jet, one_adjoint, pair_adjoint, generator,
      zero_one, zero_pair, zero_one);
  require_close(framed_adjoint_jet.overlap, raw_adjoint_jet.overlap, 2e-8,
      "framed overlap-adjoint direction differs");
  require_close(framed_adjoint_jet.one_electron, raw_adjoint_jet.one_electron, 2e-8,
      "framed one-electron-adjoint direction differs");
  require_close(xmvb::integrals::pullback_pair_kernel(
      framed_adjoint_jet.pair_kernel, q), raw_adjoint_jet.pair_kernel, 2e-8,
      "framed pair-adjoint direction differs");
}

}  // namespace

int main() {
  try {
    check_frame(false);
    check_frame(true);
    std::cout << "fused exact/RI integral frame: passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
