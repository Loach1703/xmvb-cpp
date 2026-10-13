#pragma once

#include <vector>
#include <Eigen/Core>

namespace xmvb::vb {
struct VbScfInput;

/** @brief Only the mean-field density diagonals consumed by the orbital model. */
struct DensityDiagonals {
  Eigen::MatrixXd core, active;
  /** @brief Fused core Coulomb projection for directed CSR integrals. */
  Eigen::MatrixXd core_coulomb;
};

/** @brief Exact source-span contraction without spectral density truncation. */
DensityDiagonals contract_diagonal_densities(
    const VbScfInput& input, const Eigen::MatrixXd& frame, int n_core,
    const Eigen::MatrixXd& dm1, const std::vector<double>& dm2,
    const std::vector<double>& active_coulomb);
}  // namespace xmvb::vb
