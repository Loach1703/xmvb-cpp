#include "core/workspace.hpp"
#include "dft/grid/atomic.hpp"
#include "dft/integration/potential.hpp"
#include "dft/xc/functional.hpp"
#include "integrals/ao/libcint_input.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

extern "C" {
#include "cint.h"
}

namespace {

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

template <class Callable>
void require_invalid(Callable&& call, const char* message) {
  try {
    call();
  } catch (const std::invalid_argument&) {
    return;
  }
  throw std::runtime_error(message);
}

struct RestoreWorkspace {
  std::size_t original = xmvb::workspace::limit();
  ~RestoreWorkspace() { xmvb::workspace::set_limit(original); }
};

xmvb::integrals::LibcintInput make_two_s_basis() {
  xmvb::integrals::LibcintInput input;
  input.n_atoms = 1;
  input.n_shells = 2;
  input.n_gaussian_primitives = 2;
  input.atm.assign(ATM_SLOTS, 0);
  input.atm[CHARGE_OF] = 1;
  input.atm[PTR_COORD] = PTR_ENV_START;
  input.env.assign(PTR_ENV_START + 3, 0.0);
  input.bas.assign(2 * BAS_SLOTS, 0);
  input.basidx = {0, 1, 1, 1};
  for (int shell = 0; shell < 2; ++shell) {
    const int offset = static_cast<int>(input.env.size());
    input.env.push_back(shell == 0 ? 0.7 : 1.4);
    input.env.push_back(1.0);
    input.bas[shell * BAS_SLOTS + ATOM_OF] = 0;
    input.bas[shell * BAS_SLOTS + ANG_OF] = 0;
    input.bas[shell * BAS_SLOTS + NPRIM_OF] = 1;
    input.bas[shell * BAS_SLOTS + NCTR_OF] = 1;
    input.bas[shell * BAS_SLOTS + PTR_EXP] = offset;
    input.bas[shell * BAS_SLOTS + PTR_COEFF] = offset + 1;
  }
  return input;
}

xmvb::dft::Grid make_grid() {
  xmvb::dft::Grid grid;
  grid.points.resize(12, 3);
  grid.weights.resize(12);
  for (Eigen::Index i = 0; i < 12; ++i) {
    grid.points(i, 0) = 0.13 * static_cast<double>(i - 5);
    grid.points(i, 1) = 0.07 * static_cast<double>((i * 3) % 7 - 3);
    grid.points(i, 2) = 0.11 * static_cast<double>((i * 5) % 9 - 4);
    grid.weights[i] = 0.04 + 0.002 * static_cast<double>(i);
  }
  return grid;
}

void check_dirac_exchange() {
  xmvb::dft::LdaFunctional exchange("LDA_X");
  Eigen::VectorXd rho(3);
  rho << 0.01, 0.5, 2.0;
  Eigen::VectorXd exc, vrho;
  exchange.evaluate(rho, &exc, &vrho);
  const double coefficient = std::cbrt(3.0 / std::acos(-1.0));
  for (Eigen::Index i = 0; i < rho.size(); ++i) {
    const double scale = coefficient * std::cbrt(rho[i]);
    require(std::abs(exc[i] + 0.75 * scale) < 2.0e-12,
            "Dirac exchange energy per particle mismatch");
    require(std::abs(vrho[i] + scale) < 2.0e-12,
            "Dirac exchange potential mismatch");
  }
  const Eigen::VectorXd exchange_exc = exc;
  const Eigen::VectorXd exchange_vrho = vrho;
  xmvb::dft::LdaFunctional correlation("LDA_C_VWN");
  Eigen::VectorXd correlation_exc, correlation_vrho;
  correlation.evaluate(rho, &correlation_exc, &correlation_vrho);
  xmvb::dft::LdaFunctional default_functional;
  default_functional.evaluate(rho, &exc, &vrho);
  require((exc - exchange_exc - correlation_exc).norm() < 1.0e-13 &&
              (vrho - exchange_vrho - correlation_vrho).norm() < 1.0e-13,
          "default LDA composition mismatch");

  require_invalid([] { xmvb::dft::LdaFunctional bad("GGA_X_PBE"); },
                  "GGA accepted as LDA");
  require_invalid([] { xmvb::dft::LdaFunctional bad("LDA_K_TF"); },
                  "kinetic functional accepted as XC");
  require_invalid([] { xmvb::dft::LdaFunctional bad("LDA_X,"); },
                  "empty XC component accepted");
  require_invalid([] { xmvb::dft::LdaFunctional bad("NOT_A_LIBXC_NAME"); },
                  "unknown XC functional accepted");
  rho[0] = -1.0e-5;
  require_invalid([&] { exchange.evaluate(rho, &exc, &vrho); },
                  "negative density accepted");
  rho[0] = std::numeric_limits<double>::quiet_NaN();
  require_invalid([&] { exchange.evaluate(rho, &exc, &vrho); },
                  "nonfinite density accepted");
}

void check_grid_derivative_and_tiling() {
  RestoreWorkspace restore;
  const auto input = make_two_s_basis();
  const auto grid = make_grid();
  const xmvb::dft::LdaFunctional exchange("LDA_X");
  const xmvb::dft::LdaIntegrator integrator(input, grid, exchange);
  Eigen::Matrix2d density;
  density << 0.9, 0.1, 0.1, 0.5;
  Eigen::Matrix2d direction;
  direction << 0.02, -0.01, -0.01, -0.015;

  xmvb::workspace::set_limit(4096);
  const auto full = integrator.evaluate(density);
  const double step = 1.0e-5;
  const double plus = integrator.evaluate(density + step * direction).energy;
  const double minus = integrator.evaluate(density - step * direction).energy;
  const double finite_difference = (plus - minus) / (2.0 * step);
  const double contraction =
      (full.potential.array() * direction.array()).sum();
  require(std::abs(finite_difference - contraction) <
              2.0e-7 * std::max(1.0, std::abs(contraction)),
          "XC potential is not the energy derivative");
  require(full.electron_count > 0.0 && full.energy < 0.0,
          "unexpected grid density or exchange energy");

  // nbf=2: output=32 bytes and one grid row needs (2*2+6)*8=80 bytes.
  xmvb::workspace::set_limit(192);
  const auto tiled = integrator.evaluate(density);
  require(std::abs(full.energy - tiled.energy) < 1.0e-13,
          "XC energy depends on tile size");
  require(std::abs(full.electron_count - tiled.electron_count) < 1.0e-13,
          "electron count depends on tile size");
  require((full.potential - tiled.potential).norm() < 1.0e-12,
          "XC potential depends on tile size");

  xmvb::workspace::set_limit(111);
  bool exhausted = false;
  try {
    (void)integrator.evaluate(density);
  } catch (const std::runtime_error& error) {
    exhausted = std::string(error.what()).find("resource_limit") !=
                std::string::npos;
  }
  require(exhausted, "insufficient XC workspace was not reported");
}

}  // namespace

int main() {
  check_dirac_exchange();
  check_grid_derivative_and_tiling();
  std::cout << "dft_xc=PASS\n";
}
