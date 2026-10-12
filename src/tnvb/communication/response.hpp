#pragma once

#include <vector>
#include <Eigen/Core>

namespace xmvb::tnvb {

/** Owned column-major numerical response; no determinant cache or rank-6 data. */
struct WorkerResponse {
  Eigen::MatrixXd hamiltonian;
  Eigen::MatrixXd overlap;
  Eigen::VectorXd energies;
  Eigen::MatrixXd coefficients;
  Eigen::MatrixXd density_one;
  std::vector<double> density_two_buffer;
  Eigen::MatrixXd contracted_three;
  double fit_energy_error = 0.0;
  double rdm_identity_error = 0.0;
  double hs_wall_seconds = 0.0;
  double rdm_wall_seconds = 0.0;
};

}  // namespace xmvb::tnvb
