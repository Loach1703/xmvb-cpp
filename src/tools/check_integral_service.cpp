#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "input/deck/model.hpp"
#include "input/deck/primary_basis.hpp"
#include "libcint/materialized_input.hpp"
#include "libcint/materialized_provider.hpp"
#include "parallel/mpi.hpp"
#include "guess/closed_shell_fock.hpp"
#include "vbscf/integrals/active/two_electron/construction/builder.hpp"
#include "vbscf/integrals/active/two_electron/response/adjoint.hpp"
#include "vbscf/integrals/active/two_electron/transformation/ao_pair_operator.hpp"

namespace {

/** @brief Compare every operator component, not just a scalar bilinear probe. */
void check(const double* actual, const double* reference, std::size_t count) {
  double error = 0, scale = 1;
  for (std::size_t i = 0; i < count; ++i) {
    if (!std::isfinite(actual[i]) || !std::isfinite(reference[i]))
      throw std::runtime_error("nonfinite MPI integral action");
    error = std::max(error, std::abs(actual[i] - reference[i]));
    scale = std::max(scale, std::abs(reference[i]));
  }
  std::cout << "maximum relative component error: " << error / scale << '\n';
  if (error > 2e-12 * scale)
    throw std::runtime_error("MPI integral action differs from complete reference");
}

/** @brief Independent eightfold-symmetry 2J-K reference, without the H1E stencil. */
Eigen::MatrixXd fock_reference(const xmvb::vb::AoIntegralInput& ao,
                             const Eigen::MatrixXd& density) {
  Eigen::MatrixXd f = ao.ao_core_hamiltonian_matrix;
  const auto& graph = ao.pair_graph;
  graph.for_each_integral([&](double value, int i, int j, int k, int l) {
    std::array<std::array<int, 4>, 8> permutations{{
        {i,j,k,l}, {j,i,k,l}, {i,j,l,k}, {j,i,l,k},
        {k,l,i,j}, {l,k,i,j}, {k,l,j,i}, {l,k,j,i}}};
    std::sort(permutations.begin(), permutations.end());
    const auto end = std::unique(permutations.begin(), permutations.end());
    for (auto p = permutations.begin(); p != end; ++p) {
      const auto [a,b,c,d] = *p;
      f(a,b) += 2 * value * density(c,d);
      f(a,c) -= value * density(b,d);
    }
  });
  return (0.5 * (f + f.transpose())).eval();
}

}  // namespace

/** @brief Small-fixture transport regression; only this test builds a full root reference. */
int main(int argc, char** argv) {
  try {
    xmvb::parallel::MpiSession session(argc, argv);
    if (session.run_workers()) return 0;
    if (argc != 2) throw std::invalid_argument("usage: check_integral_service input.xmi");
    const auto deck = xmvb::vb::parse_input_deck_model(argv[1]);
    const auto input = xmvb::vb::build_input_deck_basis(deck, deck.metadata.basis_name).libcint_input;
    // Limit the independent full-integral reference to a small diagnostic fixture.
    if (input.basidx.size() < 2 ||
        input.basidx[input.basidx.size() - 2] + input.basidx.back() > 64)
      throw std::invalid_argument("MPI service regression requires at most 64 AOs");
    auto distributed = session.build(input);
    auto reference = xmvb::vb::build_materialized_ao_integral_input(
        xmvb::vb::LibcintMaterializedIntegralProvider().build(input));
    const int n_bf = reference.n_basis_functions;
    const int n_pairs = n_bf * (n_bf + 1) / 2;
    const int rows = n_bf * n_bf;
    std::vector<double> source(rows), adjoint(rows), forward, transpose;
    for (int i = 0; i < rows; ++i) {
      source[i] = std::sin((i + 1) * 0.13);
      adjoint[i] = std::cos((i + 1) * 0.17);
    }
    auto expected = xmvb::vb::apply_ao_h1e(source.data(), reference, 1);
    auto actual = session.h1e(distributed, source.data(), false);
    check(actual.data(), expected.data(), rows);
    expected = xmvb::vb::apply_ao_h1e_transpose(adjoint.data(), reference, 1);
    actual = session.h1e(distributed, adjoint.data(), true);
    check(actual.data(), expected.data(), rows);
    session.fused(distributed, source.data(), adjoint.data(), &forward, &transpose);
    check(transpose.data(), expected.data(), rows);
    expected = xmvb::vb::apply_ao_h1e(source.data(), reference, 1);
    check(forward.data(), expected.data(), rows);
    // Cross the transport's 8 MiB column-tile boundary, including the tail tile.
    const int columns = (8 * 1024 * 1024 / sizeof(double)) / n_pairs + 3;
    xmvb::vb::PairMatrix coefficients(n_pairs, columns), image, full_image;
    for (Eigen::Index i = 0; i < coefficients.size(); ++i)
      coefficients.data()[i] = std::sin((i + 1) * 0.0017);
    xmvb::vb::detail::apply_exact_ao_pair_kernel(
        reference, coefficients, n_bf, columns, &full_image);
    xmvb::vb::detail::apply_exact_ao_pair_kernel(
        distributed, coefficients, n_bf, columns, &image);
    check(image.data(), full_image.data(), image.size());
    coefficients.resize(0, 0);
    image.resize(0, 0);
    full_image.resize(0, 0);
    const Eigen::MatrixXd density = Eigen::Map<const Eigen::MatrixXd>(source.data(), n_bf, n_bf);
    const auto expected_fock = fock_reference(reference, density);
    const auto native_fock = xmvb::vb::ClosedShellFockBuilder().build(density, reference);
    const auto mpi_fock = xmvb::vb::ClosedShellFockBuilder().build(density, distributed);
    check(native_fock.data(), expected_fock.data(), rows);
    check(mpi_fock.data(), expected_fock.data(), rows);

    // Many active orbitals force both the sparse transform and nonresident
    // accepted/backward/HVP path; random dense factors do not hide cross terms.
    const int n_active = n_bf;
    Eigen::MatrixXd c(n_bf, n_active), d(n_bf, n_active);
    for (Eigen::Index i = 0; i < c.size(); ++i) {
      c.data()[i] = std::sin((i + 1) * 0.037) / n_bf;
      d.data()[i] = std::cos((i + 1) * 0.029) / n_bf;
    }
    xmvb::vb::OrbitalPreparationResult prepared;
    prepared.active_sparse_row_offsets.push_back(0);
    for (int i = 0; i < n_bf; ++i) {
      for (int j = 0; j < n_active; ++j) {
        prepared.active_sparse_orbital_indices.push_back(j);
        prepared.active_sparse_values.push_back(c(i,j));
      }
      prepared.active_sparse_row_offsets.push_back(prepared.active_sparse_values.size());
    }
    auto native_active = xmvb::vb::ActiveSpaceTwoElectronBuilder().build(reference, prepared, n_active);
    auto mpi_active = xmvb::vb::ActiveSpaceTwoElectronBuilder().build(distributed, prepared, n_active);
    if (native_active.dense_ao_pair_products.size() || mpi_active.dense_ao_pair_products.size())
      throw std::runtime_error("regression failed to exercise nonresident active transform");
    xmvb::vb::detail::build_ao_pair_coefficient_rows(reference.pair_graph, c, 0, n_pairs, &coefficients);
    xmvb::vb::detail::apply_exact_ao_pair_kernel(reference, coefficients, n_bf,
        coefficients.cols(), &full_image);
    const Eigen::MatrixXd dense_tensor = coefficients.transpose() * full_image;
    std::vector<double> packed_reference;
    for (Eigen::Index i = 0; i < dense_tensor.rows(); ++i)
      for (Eigen::Index j = 0; j <= i; ++j) packed_reference.push_back(dense_tensor(i,j));
    check(native_active.packed_active_two_electron_integrals.data(),
        packed_reference.data(), packed_reference.size());
    check(mpi_active.packed_active_two_electron_integrals.data(),
        packed_reference.data(), packed_reference.size());
    std::vector<double> weights(native_active.packed_active_two_electron_integrals.size());
    for (std::size_t i = 0; i < weights.size(); ++i) weights[i] = std::sin(i * 0.11) * 0.001;
    const auto native_cache = xmvb::vb::build_exact_packed_active_two_electron_adjoint_cache(
        weights, native_active.dense_active_coefficients, reference, n_active, native_active);
    const auto mpi_cache = xmvb::vb::build_exact_packed_active_two_electron_adjoint_cache(
        weights, mpi_active.dense_active_coefficients, distributed, n_active, mpi_active);
    const auto native_gradient = xmvb::vb::backpropagate_exact_packed_active_two_electron_gradient(weights, native_cache);
    const auto mpi_gradient = xmvb::vb::backpropagate_exact_packed_active_two_electron_gradient(weights, mpi_cache);
    check(mpi_gradient.data(), native_gradient.data(), mpi_gradient.size());
    const auto native_hvp = xmvb::vb::apply_exact_packed_active_two_electron_adjoint_hessian_vector(native_cache, d, reference);
    const auto mpi_hvp = xmvb::vb::apply_exact_packed_active_two_electron_adjoint_hessian_vector(mpi_cache, d, distributed);
    check(mpi_hvp.data(), native_hvp.data(), mpi_hvp.size());

    // Large generated images cross a row-tile boundary; workers never see
    // this complete test-only image or a complete pair coefficient matrix.
    Eigen::MatrixXd wide_c(n_bf, 80), wide_d(n_bf, 80);
    for (Eigen::Index i = 0; i < wide_c.size(); ++i) {
      wide_c.data()[i] = std::sin(i * 0.019) / n_bf;
      wide_d.data()[i] = std::cos(i * 0.023) / n_bf;
    }
    xmvb::vb::PairMatrix delta, full_delta;
    xmvb::vb::detail::apply_generated_ao_pair_rows(reference, wide_c, &wide_d,
        1, n_pairs - 2, &full_image, &full_delta);
    xmvb::vb::detail::apply_generated_ao_pair_rows(distributed, wide_c, &wide_d,
        1, n_pairs - 2, &image, &delta);
    check(image.data(), full_image.data(), image.size());
    check(delta.data(), full_delta.data(), delta.size());
    std::cout << "MPI integral service regression passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
