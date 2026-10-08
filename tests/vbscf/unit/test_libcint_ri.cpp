#include <cmath>
#include <iostream>
#include <stdexcept>

extern "C" {
#include "cint.h"
}

#include "libcint/ri_provider.hpp"
#include "libcint/shell_normalization.hpp"

namespace {

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

/** @brief One normalized s Gaussian; duplicate shells make the metric singular. */
xmvb::vb::LibcintInput make_input(int n_shells) {
  xmvb::vb::LibcintInput input;
  input.n_atoms = 1;
  input.n_shells = n_shells;
  input.n_gaussian_primitives = n_shells;
  input.atm.assign(ATM_SLOTS, 0);
  input.atm[CHARGE_OF] = 1;
  input.atm[PTR_COORD] = PTR_ENV_START;
  input.env.assign(PTR_ENV_START + 3, 0.0);
  input.bas.assign(n_shells * BAS_SLOTS, 0);
  input.basidx.resize(2 * n_shells);
  for (int shell = 0; shell < n_shells; ++shell) {
    const int offset = static_cast<int>(input.env.size());
    input.env.push_back(1.0);
    input.env.push_back(xmvb::vb::shell_normalization(0, 1.0));
    input.bas[shell * BAS_SLOTS + NPRIM_OF] = 1;
    input.bas[shell * BAS_SLOTS + NCTR_OF] = 1;
    input.bas[shell * BAS_SLOTS + PTR_EXP] = offset;
    input.bas[shell * BAS_SLOTS + PTR_COEFF] = offset + 1;
    input.basidx[2 * shell] = shell;
    input.basidx[2 * shell + 1] = 1;
  }
  return input;
}

}  // namespace

int main() {
  try {
    constexpr double pi = 3.14159265358979323846;
    using xmvb::vb::radial_cartesian_scale;
    require(radial_cartesian_scale(0,0,0) == 1.0, "s normalization");
    require(radial_cartesian_scale(1,0,0) == 1.0, "p normalization");
    require(std::abs(std::pow(radial_cartesian_scale(2,0,0),2) - 4*pi/5) < 1e-14,
            "Cartesian dxx radial norm");
    require(std::abs(std::pow(radial_cartesian_scale(1,1,0),2) - 4*pi/15) < 1e-14,
            "Cartesian dxy radial norm");
    require(std::abs(std::pow(radial_cartesian_scale(3,0,0),2) - 4*pi/7) < 1e-14,
            "Cartesian fxxx radial norm");

    const auto primary = make_input(1);
    xmvb::vb::LibcintRiIntegralProvider provider;
    const auto positive = provider.build(primary, make_input(1), {});
    const auto singular = provider.build(primary, make_input(3), {});
    require(positive.n_auxiliary_functions == 1, "positive metric rank");
    require(singular.n_auxiliary_functions == 1, "duplicate auxiliary shells must be removed");
    require(singular.auxiliary_metric_matrix.rows() == 3, "original metric size");
    require(singular.metric_whitened_ao_pair_factors.rows() == 1, "rectangular factor rank");
    const double a = positive.metric_whitened_ao_pair_factors.squaredNorm();
    const double b = singular.metric_whitened_ao_pair_factors.squaredNorm();
    require(std::abs(a-b) < 1e-12, "Cholesky and spectral RI contractions differ");
    std::cout << "RI metric normalization and rank regression passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
