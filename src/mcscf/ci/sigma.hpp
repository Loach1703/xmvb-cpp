#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "mcscf/ci/integrals.hpp"

namespace xmvb::mcscf {

/** @brief Integral adjoint of one orthonormal direct-CI bilinear form. */
struct DirectCiIntegralAdjoint {
  Eigen::MatrixXd one_electron;
  Eigen::MatrixXd pair_kernel;
};

/** @brief Sigma and integral adjoint sharing one spin-string tile sweep. */
struct DirectCiResponse {
  Eigen::MatrixXd sigma;
  DirectCiIntegralAdjoint integrals;
};

/** @brief One generator action; inputs are views, output is caller-owned. */
struct DirectCiGeneratorRequest {
  Eigen::Ref<const Eigen::MatrixXd> coefficients;
  Eigen::Ref<const Eigen::MatrixXd> generator;
  Eigen::MatrixXd* result;
};

/**
 * @brief Exact Slater--Condon sigma action in an orthonormal active basis.
 */
class DirectCiSigmaAction {
public:
  DirectCiSigmaAction(
      const std::vector<std::vector<int>>& alpha_determinants,
      const std::vector<std::vector<int>>& beta_determinants,
      const ActiveIntegrals& integrals);

  /** @brief Reuses immutable spin-string links for another integral kernel.
   *
   * Orbitals must have the same ordering and count as the reference action.
   * Only the numerical integral kernels are rebuilt; links are shared, not copied.
   */
  DirectCiSigmaAction(
      const DirectCiSigmaAction& reference,
      const ActiveIntegrals& integrals);

  /**
   * @brief Applies the Hamiltonian to packed determinant-product blocks.
   *
   * The input shape is `(n_alpha, block_width * n_beta)`; beta columns for one
   * block vector are contiguous.
   */
  Eigen::MatrixXd apply(
      const Eigen::Ref<const Eigen::MatrixXd>& coefficients) const;

  /** @brief Stream the exact Hamiltonian action through beta-column callbacks.
   *
   * Each callback receives the requested alpha row indices and a
   * `(alpha_rows.size(), block_width)` column slab. The
   * reader fills its slab from an immutable input carrier; the accumulator
   * adds its slab to a caller-owned, initially zero output carrier. Neither
   * callback may alias the two temporary slabs across calls. Shared alpha
   * excitation halos and worker scatter/reduction tiles are sized together
   * under a workspace target; workers never copy the complete carrier.
   */
  void apply_tiled(
      int block_width,
      const std::function<void(int, const std::vector<int>&, Eigen::Ref<Eigen::MatrixXd>)>& read_column,
      const std::function<void(int, const std::vector<int>&, const Eigen::Ref<const Eigen::MatrixXd>&)>& add_column) const;

  /** @brief Hamiltonian diagonal on the unique alpha/beta product carrier.
   * Uses occupation lists and the retained two-electron kernel, not sigma
   * probes. The one-electron integrals must belong to this immutable action.
   */
  Eigen::MatrixXd diagonal(const Eigen::Ref<const Eigen::MatrixXd>& one_electron) const;

  /**
   * @brief Applies a general one-body orbital generator.
   *
   * The generator follows `sum_pq kappa(p,q) a_p^+ a_q` and need not be
   * symmetric. The coefficient packing is identical to `apply()`.
   */
  Eigen::MatrixXd apply_one_body_generator(
      const Eigen::Ref<const Eigen::MatrixXd>& coefficients,
      const Eigen::Ref<const Eigen::MatrixXd>& generator) const;

  /** @brief Applies equally shaped generator requests in one parallel team.
   * Inputs stay shared; each request must have a distinct, non-aliasing output.
   */
  void apply_one_body_generators(
      std::initializer_list<DirectCiGeneratorRequest> requests) const;

  /**
   * @brief Differentiates `left.dot(H * right)` with respect to the integrals.
   *
   * Both returned matrices use the symmetric full-matrix Frobenius convention:
   * their inner products with symmetric integral directions equal the exact
   * directional derivative of the bilinear form. The traversal uses the same
   * Slater--Condon connection graph as `apply()`.
   */
  DirectCiIntegralAdjoint integral_adjoint(
      const Eigen::Ref<const Eigen::MatrixXd>& left,
      const Eigen::Ref<const Eigen::MatrixXd>& right) const;

  /** @brief Returns H * response and the integral derivative of response.dot(H * state).
   * Reuses the response scatter; two shared tiles serve both contractions.
   * Inputs follow apply() packing and may be non-owning strided blocks.
   */
  DirectCiResponse apply_with_adjoint(
      const Eigen::Ref<const Eigen::MatrixXd>& response,
      const Eigen::Ref<const Eigen::MatrixXd>& state) const;

  /**
   * @brief Unsymmetrized one-body adjoint for an orbital-generator action.
   *
   * The returned matrix `D` satisfies
   * `left.dot(E(kappa) * right) = D.cwiseProduct(kappa).sum()` for a general,
   * not necessarily symmetric one-body generator `kappa`.
   */
  Eigen::MatrixXd one_body_generator_adjoint(
      const Eigen::Ref<const Eigen::MatrixXd>& left,
      const Eigen::Ref<const Eigen::MatrixXd>& right) const;

  /** @brief Numerical matrices plus referenced links, not exclusively owned bytes. */
  std::size_t dynamic_bytes() const noexcept;

private:
  struct DensityConnections {
    std::vector<int> sources;
    std::vector<int> pairs;
    std::vector<double> signs;
    std::vector<int> created_orbitals;
    std::vector<int> annihilated_orbitals;

    std::size_t size() const noexcept { return sources.size(); }
    void reserve(std::size_t capacity);
    void append(
        int source,
        int pair,
        double sign,
        int created_orbital,
        int annihilated_orbital);
    std::size_t dynamic_bytes() const noexcept;
  };

  struct SpinConnections {
    std::vector<DensityConnections> singles;
    std::vector<std::vector<int>> occupied;
  };

  SpinConnections build_spin_connections(
      const std::vector<std::vector<int>>& determinants) const;
  Eigen::MatrixXd build_hamiltonian_pair_kernel(
      const ActiveIntegrals& integrals,
      int n_electrons) const;
  /** @brief Scatter into (block_width * alpha_count, orbital_pairs), strings contiguous. */
  void scatter_one_body_intermediate(
      const Eigen::Ref<const Eigen::MatrixXd>& coefficients,
      int alpha_begin,
      int alpha_count,
      int beta_target,
      Eigen::Ref<Eigen::MatrixXd> intermediate) const;
  const SpinConnections& beta_connections() const noexcept;

  /** @brief Shared adjoint sweep; non-null sigma additionally receives H * left. */
  DirectCiIntegralAdjoint contract_adjoint(
      const Eigen::Ref<const Eigen::MatrixXd>& left,
      const Eigen::Ref<const Eigen::MatrixXd>& right,
      Eigen::MatrixXd* sigma) const;

  /** @brief Cooperative sigma contraction inside the current OpenMP team.
   * Each output beta column has one owner. All workers must call this method.
   */
  void apply_sigma_tile(
      int alpha_begin, int alpha_count, int block_begin, int block_count,
      int beta_begin, int beta_count,
      const Eigen::Ref<const Eigen::MatrixXd>& scattered,
      Eigen::Ref<Eigen::MatrixXd> workspace,
      Eigen::Ref<Eigen::MatrixXd> sigma) const;

  int n_orbitals_ = 0;
  int n_alpha_ = 0;
  int n_beta_ = 0;
  Eigen::MatrixXd pair_kernel_;
  Eigen::MatrixXd hamiltonian_pair_kernel_;
  std::shared_ptr<const SpinConnections> alpha_;
  std::shared_ptr<const SpinConnections> distinct_beta_;
};

}  // namespace xmvb::mcscf
