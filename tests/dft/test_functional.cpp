#include "dft/xc/functional.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace {

using xmvb::dft::Functional;
using xmvb::dft::XcValues;

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance,
          const char* message) {
  if (!(std::abs(actual - expected) <= tolerance)) {
    throw std::runtime_error(std::string(message) + ": actual " +
                             std::to_string(actual) + ", expected " +
                             std::to_string(expected));
  }
}

template <class Action>
void require_invalid(Action&& action, const char* message) {
  try {
    action();
  } catch (const std::invalid_argument&) {
    return;
  }
  throw std::runtime_error(message);
}

XcValues evaluate(const Functional& functional, const Eigen::MatrixXd& rho,
                  const Eigen::MatrixXd& sigma) {
  XcValues values;
  functional.evaluate(rho, sigma, &values);
  return values;
}

double total_energy(const Functional& functional, const Eigen::MatrixXd& rho,
                    const Eigen::MatrixXd& sigma) {
  const XcValues values = evaluate(functional, rho, sigma);
  double energy = 0.0;
  for (Eigen::Index i = 0; i < rho.rows(); ++i) {
    energy += rho.row(i).sum() * values.exc[i];
  }
  return energy;
}

void check_lda_exchange() {
  static_assert(!std::is_copy_constructible_v<Functional>);
  static_assert(std::is_move_constructible_v<Functional>);
  const double coefficient = std::cbrt(3.0 / std::acos(-1.0));
  const double spin_factor = std::cbrt(2.0);

  Functional restricted("LDA_X");
  require(restricted.n_spin() == 1 && !restricted.is_gga(),
          "LDA restricted metadata");
  Eigen::MatrixXd rho(3, 1);
  rho << 0.01, 0.5, 2.0;
  const Eigen::MatrixXd no_sigma(3, 0);
  const auto values = evaluate(restricted, rho, no_sigma);
  require(values.exc.size() == 3 && values.vrho.rows() == 3 &&
              values.vrho.cols() == 1 && values.vsigma.rows() == 3 &&
              values.vsigma.cols() == 0,
          "LDA output shapes");
  for (Eigen::Index i = 0; i < rho.rows(); ++i) {
    const double scale = coefficient * std::cbrt(rho(i, 0));
    near(values.exc[i], -0.75 * scale, 2e-12,
         "restricted Dirac exchange energy");
    near(values.vrho(i, 0), -scale, 2e-12,
         "restricted Dirac exchange potential");
  }

  Functional polarized("LDA_X", 2);
  require(polarized.n_spin() == 2 && !polarized.is_gga(),
          "LDA polarized metadata");
  Eigen::MatrixXd spin_rho(3, 2);
  spin_rho << 0.4, 0.2, 0.7, 0.0, 0.03, 1.2;
  const auto spin_values = evaluate(polarized, spin_rho,
                                    Eigen::MatrixXd(3, 0));
  for (Eigen::Index i = 0; i < spin_rho.rows(); ++i) {
    const double alpha = spin_rho(i, 0);
    const double beta = spin_rho(i, 1);
    const double expected_exc =
        -0.75 * spin_factor * coefficient *
        (alpha * std::cbrt(alpha) + beta * std::cbrt(beta)) /
        (alpha + beta);
    near(spin_values.exc[i], expected_exc, 2e-12,
         "polarized Dirac exchange energy");
    near(spin_values.vrho(i, 0),
         -spin_factor * coefficient * std::cbrt(alpha), 2e-12,
         "alpha Dirac exchange potential");
    near(spin_values.vrho(i, 1),
         -spin_factor * coefficient * std::cbrt(beta), 2e-12,
         "beta Dirac exchange potential");
  }
  Eigen::MatrixXd full_rho(1, 1);
  full_rho(0, 0) = 0.7;
  near(spin_values.exc[1],
       spin_factor * evaluate(restricted, full_rho, Eigen::MatrixXd(1, 0))
                         .exc[0],
       2e-12, "full-polarization exchange scaling");

  Eigen::MatrixXd ignored_sigma = Eigen::MatrixXd::Constant(3, 1, 42.0);
  const auto ignored = evaluate(restricted, rho, ignored_sigma);
  near((ignored.exc - values.exc).norm(), 0.0, 1e-14,
       "LDA ignores sigma energy");
  near((ignored.vrho - values.vrho).norm(), 0.0, 1e-14,
       "LDA ignores sigma potential");
}

void check_composition() {
  Eigen::MatrixXd rho(2, 2);
  rho << 0.7, 0.2, 0.3, 0.8;
  Eigen::MatrixXd sigma(2, 3);
  sigma << 0.04, -0.01, 0.07, 0.08, 0.01, 0.03;
  Functional lda("LDA_X", 2);
  Functional gga("GGA_C_PBE", 2);
  Functional mixed(" LDA_X , GGA_C_PBE ", 2);
  require(mixed.is_gga(), "mixed LDA/GGA metadata");
  const auto a = evaluate(lda, rho, Eigen::MatrixXd(2, 0));
  const auto b = evaluate(gga, rho, sigma);
  const auto sum = evaluate(mixed, rho, sigma);
  near((sum.exc - a.exc - b.exc).norm(), 0.0, 1e-12,
       "mixed XC energy composition");
  near((sum.vrho - a.vrho - b.vrho).norm(), 0.0, 1e-12,
       "mixed XC density derivative composition");
  near((sum.vsigma - b.vsigma).norm(), 0.0, 1e-12,
       "mixed XC sigma derivative composition");
}

void check_derivatives(const Functional& functional, Eigen::MatrixXd rho,
                       Eigen::MatrixXd sigma) {
  const auto values = evaluate(functional, rho, sigma);
  const double step = 1e-6;
  for (Eigen::Index i = 0; i < rho.rows(); ++i) {
    for (Eigen::Index channel = 0; channel < rho.cols(); ++channel) {
      rho(i, channel) += step;
      const double plus = total_energy(functional, rho, sigma);
      rho(i, channel) -= 2.0 * step;
      const double minus = total_energy(functional, rho, sigma);
      rho(i, channel) += step;
      const double derivative = (plus - minus) / (2.0 * step);
      near(derivative, values.vrho(i, channel),
           2e-7 * std::max(1.0, std::abs(derivative)),
           "Libxc density finite difference");
    }
    for (Eigen::Index channel = 0; channel < sigma.cols(); ++channel) {
      sigma(i, channel) += step;
      const double plus = total_energy(functional, rho, sigma);
      sigma(i, channel) -= 2.0 * step;
      const double minus = total_energy(functional, rho, sigma);
      sigma(i, channel) += step;
      const double derivative = (plus - minus) / (2.0 * step);
      near(derivative, values.vsigma(i, channel),
           2e-7 * std::max(1.0, std::abs(derivative)),
           "Libxc sigma finite difference");
    }
  }
}

void check_gga() {
  Functional restricted("GGA_X_PBE,GGA_C_PBE");
  require(restricted.is_gga() && restricted.n_spin() == 1,
          "restricted PBE metadata");
  Eigen::MatrixXd rho(3, 1), sigma(3, 1);
  rho << 0.3, 0.8, 1.2;
  sigma << 0.017, 0.09, 0.21;
  check_derivatives(restricted, rho, sigma);

  Functional polarized("GGA_X_PBE,GGA_C_PBE", 2);
  Eigen::MatrixXd spin_rho(3, 2), spin_sigma(3, 3);
  spin_rho << 0.4, 0.2, 0.8, 0.3, 0.35, 0.9;
  spin_sigma << 0.018, 0.005, 0.025,
                0.06, -0.01, 0.04,
                0.13, -0.02, 0.09;
  check_derivatives(polarized, spin_rho, spin_sigma);
  check_derivatives(Functional("LDA_X,GGA_C_PBE", 2), spin_rho, spin_sigma);

  const auto batch = evaluate(polarized, spin_rho, spin_sigma);
  for (Eigen::Index i = 0; i < spin_rho.rows(); ++i) {
    Eigen::MatrixXd one_rho = spin_rho.row(i);
    Eigen::MatrixXd one_sigma = spin_sigma.row(i);
    const auto one = evaluate(polarized, one_rho, one_sigma);
    near(batch.exc[i], one.exc[0], 1e-13, "interleaved XC energy");
    near((batch.vrho.row(i) - one.vrho.row(0)).norm(), 0.0, 1e-13,
         "interleaved XC vrho");
    near((batch.vsigma.row(i) - one.vsigma.row(0)).norm(), 0.0, 1e-13,
         "interleaved XC vsigma");
  }

  Eigen::MatrixXd edge_rho(2, 2), edge_sigma(2, 3);
  edge_rho << 0.7, 0.0, 0.2, 0.1;
  edge_sigma << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
  const auto edge = evaluate(polarized, edge_rho, edge_sigma);
  require(edge.exc.allFinite() && edge.vrho.allFinite() &&
              edge.vsigma.allFinite(),
          "zero spin density and gradient must remain finite");
}

void check_rejections() {
  require_invalid([] { Functional bad("LDA_X", 0); }, "n_spin=0 accepted");
  require_invalid([] { Functional bad("LDA_X", 3); }, "n_spin=3 accepted");
  require_invalid([] { Functional bad("LDA_K_TF"); },
                  "kinetic component accepted");
  require_invalid([] { Functional bad("HYB_GGA_XC_PBEH"); },
                  "hybrid component accepted");
  require_invalid([] { Functional bad("MGGA_X_TPSS"); },
                  "meta-GGA component accepted");
  require_invalid([] { Functional bad("NOT_A_LIBXC_NAME"); },
                  "unknown component accepted");
  require_invalid([] { Functional bad("LDA_X,"); },
                  "empty component accepted");

  Functional gga("GGA_X_PBE", 2);
  Eigen::MatrixXd rho = Eigen::MatrixXd::Constant(1, 2, 0.5);
  Eigen::MatrixXd sigma = Eigen::MatrixXd::Zero(1, 3);
  XcValues values;
  require_invalid([&] { gga.evaluate(rho, Eigen::MatrixXd(1, 0), &values); },
                  "missing GGA sigma accepted");
  rho(0, 0) = -1e-4;
  require_invalid([&] { gga.evaluate(rho, sigma, &values); },
                  "negative density accepted");
  rho(0, 0) = std::numeric_limits<double>::quiet_NaN();
  require_invalid([&] { gga.evaluate(rho, sigma, &values); },
                  "nonfinite density accepted");
  rho(0, 0) = 0.5;
  sigma(0, 2) = -1e-4;
  require_invalid([&] { gga.evaluate(rho, sigma, &values); },
                  "negative diagonal sigma accepted");
  sigma(0, 2) = 0.0;
  sigma(0, 1) = std::numeric_limits<double>::infinity();
  require_invalid([&] { gga.evaluate(rho, sigma, &values); },
                  "nonfinite cross sigma accepted");
  sigma(0, 0) = 0.02;
  sigma(0, 2) = 0.02;
  sigma(0, 1) = -0.01;
  require(evaluate(gga, rho, sigma).exc.allFinite(),
          "negative cross sigma rejected");

  rho(0, 1) = -1e-15;
  require(evaluate(Functional("LDA_X", 2), rho, Eigen::MatrixXd(1, 0))
              .exc.allFinite(),
          "roundoff-scale negative density rejected");
  rho(0, 1) = 0.5;
  sigma(0, 0) = -1e-16;
  sigma(0, 1) = 0.0;
  const auto rounded = evaluate(gga, rho, sigma);
  require(rounded.exc.allFinite() && rounded.vrho.allFinite() &&
              rounded.vsigma.allFinite(),
          "roundoff-scale negative sigma rejected");
}

}  // namespace

int main() {
  try {
    check_lda_exchange();
    check_composition();
    check_gga();
    check_rejections();
  } catch (const std::exception& error) {
    std::cerr << "DFT functional test failed: " << error.what() << '\n';
    return 1;
  }
  std::cout << "dft_functional=PASS\n";
}
