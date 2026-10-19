#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <stdexcept>

#include "dft/grid/atomic.hpp"
#include "dft/integration/potential.hpp"
#include "dft/xc/functional.hpp"
#include "input/deck/model.hpp"
#include "input/loading/molecule.hpp"
#include "integrals/basis/values.hpp"
#include "linalg/eigensolver.hpp"

namespace {
void write_matrix(const std::filesystem::path& path,
    const Eigen::MatrixXd& matrix) {
  std::ofstream stream(path, std::ios::binary);
  const std::int64_t shape[] = {matrix.rows(), matrix.cols()};
  stream.write(reinterpret_cast<const char*>(shape), sizeof(shape));
  stream.write(reinterpret_cast<const char*>(matrix.data()), matrix.size() * sizeof(double));
  if (!stream) throw std::runtime_error("cannot write DFT oracle matrix");
}
}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 3 && argc != 5)
      throw std::invalid_argument("expected input.xmi, oracle directory and optional radial/angular counts");
    const auto deck = xmvb::vb::parse_input_deck_model(argv[1]);
    const auto input = xmvb::input::load_molecule(deck);
    const std::filesystem::path out(argv[2]);
    std::filesystem::create_directories(out);
    const int n_bf = input.ao.n_basis_functions;
    xmvb::dft::GridOptions options;
    options.radial_points = deck.metadata.grid_radial_points;
    options.angular_points = deck.metadata.grid_angular_points;
    if (argc == 5) {
      options.radial_points = std::stoi(argv[3]);
      options.angular_points = std::stoi(argv[4]);
    }
    const auto grid = xmvb::dft::make_grid(input.libcint_input, options);
    const xmvb::integrals::BasisEvaluator basis(input.libcint_input);
    xmvb::integrals::BasisValues values;
    Eigen::MatrixXd overlap = Eigen::MatrixXd::Zero(n_bf, n_bf);
    for (Eigen::Index begin = 0; begin < grid.points.rows(); begin += 4096) {
      const Eigen::Index rows = std::min<Eigen::Index>(4096, grid.points.rows() - begin);
      basis.evaluate(grid.points.middleRows(begin, rows), 0, &values);
      const Eigen::MatrixXd weighted = grid.weights.segment(begin, rows).asDiagonal() * values.value;
      overlap.noalias() += values.value.transpose() * weighted;
    }
    const double overlap_error = (overlap - input.overlap).cwiseAbs().maxCoeff();
    std::cout << std::setprecision(17) << "overlap_error=" << overlap_error << '\n';

    xmvb::linalg::GeneralizedEigensolver solver;
    const auto& core = input.ao.ao_core_hamiltonian_matrix;
    const auto eigen = solver.solve_dense(
        std::span<const double>(core.data(), core.size()),
        std::span<const double>(input.overlap.data(), input.overlap.size()), n_bf);
    const Eigen::Map<const Eigen::MatrixXd> orbitals(eigen.eigenvector_matrix.data(), n_bf, n_bf);
    const Eigen::MatrixXd density = 2 * orbitals.leftCols(input.n_electrons / 2) *
        orbitals.leftCols(input.n_electrons / 2).transpose();
    const xmvb::dft::LdaFunctional functional(deck.metadata.xc_functional);
    const xmvb::dft::LdaIntegrator integrator(input.libcint_input, grid, functional);
    const auto xc = integrator.evaluate(density);
    Eigen::MatrixXd probes(31, 3);
    for (int row = 0; row < probes.rows(); ++row)
      for (int axis = 0; axis < 3; ++axis)
        probes(row, axis) = 0.047 * ((row * (axis + 3)) % 29 - 14);
    basis.evaluate(probes, 1, &values);
    write_matrix(out / "probes.bin", probes);
    write_matrix(out / "ao.bin", values.value);
    for (int axis = 0; axis < 3; ++axis)
      write_matrix(out / ("ao-d" + std::to_string(axis) + ".bin"), values.gradient[axis]);
    write_matrix(out / "grid.bin", grid.points);
    write_matrix(out / "weights.bin", grid.weights);
    write_matrix(out / "overlap.bin", input.overlap);
    write_matrix(out / "core.bin", core);
    write_matrix(out / "density.bin", density);
    write_matrix(out / "vxc.bin", xc.potential);
    std::ofstream summary(out / "summary.txt");
    summary << std::setprecision(17)
            << "overlap_error=" << overlap_error << '\n'
            << "xc_energy=" << xc.energy << '\n'
            << "grid_electrons=" << xc.electron_count << '\n'
            << "nuclear_energy=" << input.nuclear_energy << '\n';
    if (overlap_error > 1e-6)
      throw std::runtime_error("DFT quadrature does not reproduce AO overlap");
    std::cout << "DFT molecular quadrature and oracle export: PASS\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
