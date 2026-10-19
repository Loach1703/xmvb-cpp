#pragma once

#include <cstddef>
#include <string>
#include "dft/grid/atomic.hpp"
#include "integrals/ao/input.hpp"
#include "scf/driver/solver.hpp"

namespace xmvb::dft {

struct RksOptions {
  GridOptions grid;
  std::string functional = "LDA_X,LDA_C_VWN";
  scf::Options scf;
};

struct RksResult {
  scf::Result state;
  std::size_t grid_points = 0;
  double initial_electronic_energy = 0.0;
  double xc_energy = 0.0;
  double integrated_electrons = 0.0;
};

/** @brief Restricted LDA Kohn-Sham SCF, using the common SCF driver and exact J. */
RksResult solve_rks(int n_electrons,
    const integrals::LibcintInput& basis,
    const Eigen::Ref<const Eigen::MatrixXd>& overlap,
    const integrals::AoIntegralInput& ao,
    const RksOptions& options = {},
    const Eigen::MatrixXd& initial_density = {});

}  // namespace xmvb::dft
