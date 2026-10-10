#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

#include <Eigen/Core>

#include "vbscf/integrals/active/two_electron/construction/builder.hpp"
#include "vbscf/integrals/active/two_electron/response/adjoint.hpp"
#include "vbscf/integrals/active/two_electron/response/directional.hpp"
#include "vbscf/integrals/active/two_electron/response/internal.hpp"
#include "vbscf/integrals/active/two_electron/transformation/ao_pair_operator.hpp"
#include "vbscf/integrals/ao/pairs/two_electron_index.hpp"

namespace {

void require(bool condition, const char* message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void check_direction_block(const xmvb::vb::AoIntegralInput& integrals) {
  using namespace xmvb::vb;
  Eigen::MatrixXd orbitals(2, 2);
  orbitals << 1.0, 0.3, -0.2, 0.9;
  std::vector<Eigen::MatrixXd> directions(3, Eigen::MatrixXd::Zero(2, 2));
  directions[0] << 0.2, -0.4, 0.1, 0.3;
  directions[1] << -0.3, 0.1, 0.5, -0.2;
  ExactCtxPairMatrix pair_map;
  detail::build_ao_pair_to_active_pair_coefficients(
      orbitals, 2, 2, detail::build_active_pair_list(2), &pair_map);
  ActiveSpaceTwoElectronResult accepted;
  detail::apply_exact_ao_pair_kernel(
      integrals, pair_map, 2, 3, &accepted.dense_ao_pair_products);
  auto cache = build_exact_packed_active_two_electron_adjoint_cache(
      {0.2, -0.1, 0.3, 0.4, -0.2, 0.5}, orbitals, integrals, 2, accepted);
  std::vector<ExactCtxPairMatrix> products;
  std::vector<Eigen::MatrixXd> fixed;
  const auto batch = compute_exact_packed_active_two_electron_integral_directional_derivative_batch(
      cache, directions, integrals, &products, &fixed);
  const auto integral_only = compute_exact_packed_active_two_electron_integral_directional_derivative_batch(
      cache, directions, integrals);
  require(batch.isApprox(integral_only, 1e-13), "integral-only batch changed delta ERI");
  for (int column = 0; column < 3; ++column) {
    ExactPackedActiveTwoElectronDirectionalDerivativeWorkspace workspace;
    std::vector<double> scalar;
    compute_exact_packed_active_two_electron_integral_directional_derivative(
        cache, directions[column], integrals, &workspace, &scalar);
    require(batch.col(column).isApprox(Eigen::Map<const Eigen::VectorXd>(scalar.data(), scalar.size()), 1e-13),
            "block exact ERI differs from scalar direction");
    require(products[column].isApprox(workspace.directional_pair_products, 1e-13),
            "block half transform differs from scalar image");
    require(fixed[column].size() == 0, "resident batch returned a stale streamed adjoint");
    ExactPackedActiveTwoElectronApplyWorkspace pullback;
    Eigen::MatrixXd fused;
    apply_exact_packed_active_two_electron_adjoint_hessian_vector_fused(
        cache, directions[column], integrals, products[column], &pullback, &fused);
    require(fused.isApprox(apply_exact_packed_active_two_electron_adjoint_hessian_vector(
                cache, directions[column], integrals), 1e-13),
            "reusing block images changed orbital HVP");
  }
  if (!integrals.pair_graph.has_compact_integral_stream()) {
    cache.accepted_pair_products = nullptr;
    const auto streamed = compute_exact_packed_active_two_electron_integral_directional_derivative_batch(
        cache, directions, integrals);
    require(streamed.isApprox(batch, 1e-13), "shared accepted tiles changed delta ERI");
    const auto with_adjoint = compute_exact_packed_active_two_electron_integral_directional_derivative_batch(
        cache, directions, integrals, &products, &fixed);
    require(with_adjoint.isApprox(batch, 1e-13), "streamed orbital batch changed delta ERI");
    for (int column = 0; column < 3; ++column) {
      require(products[column].size() == 0, "streamed batch retained a full directional image");
      require(fixed[column].isApprox(apply_exact_packed_active_two_electron_adjoint_hessian_vector(
                  cache, directions[column], integrals), 1e-13),
              "streamed batch returned an incorrect fixed adjoint");
    }
  }
  const auto empty = compute_exact_packed_active_two_electron_integral_directional_derivative_batch(
      cache, {}, integrals, &products, &fixed);
  require(empty.cols() == 0 && products.empty() && fixed.empty(), "empty batch retained stale outputs");
  bool rejected = false;
  try {
    compute_exact_packed_active_two_electron_integral_directional_derivative_batch(
        cache, {Eigen::MatrixXd::Zero(1, 2)}, integrals);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  require(rejected, "batch accepted a malformed direction");
}

}  // namespace

int main() {
  constexpr int n_active = 11;

  const xmvb::vb::AoPairGraph ordered_graph =
      xmvb::vb::build_ao_pair_graph(
          {2, 1, 2}, {0, 1, 2}, {1.25, -0.5, 3.0}, n_active);
  constexpr std::array<int, 3> expected_left = {2, 1, 2};
  constexpr std::array<int, 3> expected_right = {0, 1, 2};
  constexpr std::array<double, 3> expected_values = {1.25, -0.5, 3.0};
  require(ordered_graph.integral_count() == 3, "AO-pair graph lost an integral");
  for (std::size_t index = 0; index < ordered_graph.integral_count(); ++index) {
    const xmvb::vb::AoPairIntegral integral = ordered_graph.integral(index);
    require(
        integral.left_pair == expected_left[index] &&
            integral.right_pair == expected_right[index] &&
            integral.value == expected_values[index],
        "AO-pair graph changed unique-integral order");
  }

  // The large-AO representation keeps the canonical unique stream and omits
  // the duplicated symmetric CSR edges.  Its block action must remain
  // numerically identical to the fast CSR representation.
  constexpr int compact_n_bf = 2;
  const xmvb::vb::AoPairGraph csr_graph =
      xmvb::vb::build_ao_pair_graph(
          {0, 1, 2}, {0, 0, 1}, {2.0, -0.25, 0.75}, compact_n_bf);
  xmvb::vb::AoPairGraph compact_graph;
  compact_graph.integral_rows = {0, 1, 2};
  compact_graph.integral_columns = {0, 0, 1};
  compact_graph.integral_values = {2.0, -0.25, 0.75};
  compact_graph.pair_first = csr_graph.pair_first;
  compact_graph.pair_second = csr_graph.pair_second;

  xmvb::vb::ExactCtxPairMatrix coefficients(3, 2);
  coefficients << 0.5, -1.0,
                  2.0, 0.25,
                 -0.5, 1.5;
  xmvb::vb::AoIntegralInput csr_input;
  csr_input.n_basis_functions = compact_n_bf;
  csr_input.pair_graph = csr_graph;
  xmvb::vb::AoIntegralInput compact_input;
  compact_input.n_basis_functions = compact_n_bf;
  compact_input.pair_graph = compact_graph;
  xmvb::vb::ExactCtxPairMatrix csr_products;
  xmvb::vb::ExactCtxPairMatrix compact_products;
  xmvb::vb::detail::apply_exact_ao_pair_kernel(
      csr_input, coefficients, compact_n_bf, 2, &csr_products);
  xmvb::vb::detail::apply_exact_ao_pair_kernel(
      compact_input, coefficients, compact_n_bf, 2, &compact_products);
  require(
      compact_products.isApprox(csr_products, 1.0e-14),
      "compact unique-integral action differs from symmetric CSR action");
  // Exercise SIMD lane tails and overwritten output independently of the
  // integral-transform code, including widths larger than a vector register.
  Eigen::Matrix3d reference_operator = Eigen::Matrix3d::Zero();
  for (int row = 0; row < 3; ++row)
    for (int edge = csr_graph.row_offsets[row]; edge < csr_graph.row_offsets[row + 1]; ++edge)
      reference_operator(row, csr_graph.columns[edge]) += csr_graph.values[edge];
  for (const int width : {0, 1, 2, 3, 4, 7, 8, 15, 16, 31, 32, 36, 64, 65}) {
    xmvb::vb::ExactCtxPairMatrix inputs(3, width);
    for (int row = 0; row < 3; ++row)
      for (int column = 0; column < width; ++column)
        inputs(row, column) = std::sin((row + 0.3) * (column + 0.7));
    xmvb::vb::ExactCtxPairMatrix result =
        xmvb::vb::ExactCtxPairMatrix::Constant(3, width, 1.0e20);
    xmvb::vb::detail::apply_exact_ao_pair_kernel(
        csr_input, inputs, compact_n_bf, width, &result);
    require(result.isApprox(reference_operator * inputs, 2.0e-13),
            "SIMD CSR action lost lane tails or output initialization");
  }
  check_direction_block(csr_input);
  check_direction_block(compact_input);

  xmvb::vb::OrbitalPreparationResult orbitals;
  orbitals.active_sparse_row_offsets.resize(n_active + 1);
  orbitals.active_sparse_orbital_indices.resize(n_active);
  orbitals.active_sparse_values.assign(n_active, 1.0);
  for (int index = 0; index < n_active; ++index) {
    orbitals.active_sparse_row_offsets[index] = index;
    orbitals.active_sparse_orbital_indices[index] = index;
  }
  orbitals.active_sparse_row_offsets[n_active] = n_active;

  xmvb::vb::AoIntegralInput integrals;
  integrals.n_basis_functions = n_active;
  integrals.pair_graph =
      xmvb::vb::build_ao_pair_graph({0}, {0}, {2.0}, n_active);

  const xmvb::vb::ActiveSpaceTwoElectronBuilder builder;
  const auto result = builder.build(integrals, orbitals, n_active);

  require(
      result.dense_active_coefficients.rows() == n_active &&
          result.dense_active_coefficients.cols() == n_active,
      "sparse active-ERI path omitted dense active coefficients");
  require(
      result.dense_active_coefficients.isApprox(
          Eigen::MatrixXd::Identity(n_active, n_active), 0.0),
      "sparse active-ERI path changed active coefficients");
  require(
      !result.packed_active_two_electron_integrals.empty() &&
          std::abs(result.packed_active_two_electron_integrals.front() - 2.0) <
              1.0e-14,
      "sparse active-ERI contraction returned the wrong integral");
  require(
      result.dense_ao_pair_products.size() == 0,
      "memory-bounded active-ERI path retained dense AO-pair products");
  // One off-diagonal AO edge contributes through both orientations to a
  // diagonal active-pair element: 2 * K(0,1) * B(0,0) * B(1,0) = 8.
  integrals.pair_graph = xmvb::vb::build_ao_pair_graph({0}, {1}, {2.0}, n_active);
  orbitals.active_sparse_orbital_indices[1] = 0;
  const auto mixed = builder.build(integrals, orbitals, n_active);
  require(mixed.dense_ao_pair_products.size() == 0,
          "mixed-pair regression did not exercise sparse construction");
  if (std::abs(mixed.packed_active_two_electron_integrals.front() - 8.0) >= 1e-14)
    throw std::runtime_error("sparse active diagonal: expected 8, got " +
        std::to_string(mixed.packed_active_two_electron_integrals.front()));
  // A 64-orbital packed tensor is about 16.5 MiB. Its sparse accumulator
  // must not become 32 full copies just because OMP_NUM_THREADS=32.
  constexpr int large_active = 64;
  xmvb::vb::OrbitalPreparationResult large_orbitals;
  for (int i = 0; i < large_active; ++i) {
    large_orbitals.active_sparse_row_offsets.push_back(i);
    large_orbitals.active_sparse_orbital_indices.push_back(i == 1 ? 0 : i);
    large_orbitals.active_sparse_values.push_back(1.0);
  }
  large_orbitals.active_sparse_row_offsets.push_back(large_active);
  xmvb::vb::AoIntegralInput large_integrals;
  large_integrals.n_basis_functions = large_active;
  large_integrals.pair_graph = xmvb::vb::build_ao_pair_graph({0}, {1}, {2.0}, large_active);
  const auto large_result = builder.build(large_integrals, large_orbitals, large_active);
  require(large_result.dense_ao_pair_products.size() == 0 &&
          std::abs(large_result.packed_active_two_electron_integrals.front() - 8.0) < 1e-14,
          "bounded large sparse accumulator changed the active integral");
  return 0;
}
