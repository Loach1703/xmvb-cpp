#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "input/deck/model.hpp"
#include "input/deck/primary_basis.hpp"
#include "libcint/materialized_input.hpp"
#include "libcint/materialized_provider.hpp"
#include "parallel/mpi.hpp"
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
  if (error > 2e-12 * scale)
    throw std::runtime_error("MPI integral action differs from complete reference");
  std::cout << "maximum relative component error: " << error / scale << '\n';
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
    std::cout << "MPI integral service regression passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
