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

class ProductSupport;

/** @brief Add one adjoint column range directly to a caller-owned spin column.
 * Arguments are local beta, RHS, first/last indices into alpha_rows, and a
 * destination indexed by the original alpha string, and an optional exact
 * support mask in that same indexing. A null mask requests all rows.
 * Disjoint ranges may run concurrently; the callback must complete before
 * the batch is released.
 */
using DirectCiColumnAccumulator = std::function<void(int, int, int, int, double*, const std::uint8_t*)>;

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
  /** @brief E_(spin,P) applied to single-spin source columns.
   * Output column P*width+j contains one channel image of source j.
   * No alpha-by-beta product is constructed.
   */
  void scatter_spin(bool beta,
      const Eigen::Ref<const Eigen::MatrixXd>& factors, int first, int count,
      Eigen::Ref<Eigen::MatrixXd> channels) const;

  /** @brief Accumulates sum_P E_(spin,P)^T channels[P] into single-spin images. */
  void gather_spin(bool beta, const Eigen::Ref<const Eigen::MatrixXd>& channels,
      int first, int count, Eigen::Ref<Eigen::MatrixXd> images) const;

  /** @brief Symmetric K in H = (E_alpha+E_beta)^T K (E_alpha+E_beta).
   * Includes the absorbed one-electron Hamiltonian. Borrowed immutable data.
   */
  const Eigen::MatrixXd& channel_kernel() const noexcept { return hamiltonian_pair_kernel_; }
  int n_alpha_strings() const noexcept { return n_alpha_; }
  int n_beta_strings() const noexcept { return n_beta_; }

  /** @brief Pulls a K-channel cotangent back to symmetric h/(pq|rs) adjoints.
   * Shared by product and source-factor contractions; includes the absorbed
   * one-electron and normal-ordering terms, with full Frobenius conventions.
   */
  DirectCiIntegralAdjoint pullback_channel_kernel(
      const Eigen::Ref<const Eigen::MatrixXd>& channel_adjoint) const;

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
   * block vector are contiguous. A nonzero exchange_parity requires certified
   * common transpose parity in identical spin spaces and returns a full image.
   */
  Eigen::MatrixXd apply(
      const Eigen::Ref<const Eigen::MatrixXd>& coefficients, int exchange_parity = 0) const;

  /** @brief Stream the exact Hamiltonian action through beta-column callbacks.
   *
   * Each callback receives sorted, distinct beta and alpha indices and a
   * `(alpha_rows.size(), betas.size()*block_width)` slab, with column
   * `beta_local*block_width+rhs`. The
   * reader fills its slab from an immutable input carrier. The writer receives
   * a column accumulator, which gathers directly into caller-owned output;
   * no intermediate image slab or copy is needed. Neither callback may retain
   * batch views or the accumulator across calls. Generated targets
   * span multiple beta strings; their excitation neighborhoods are merged
   * within the batch. Callbacks run in an OpenMP single region and may use
   * taskgroups, but must finish their tasks before returning. The team shares
   * bounded halos and channel slabs, not carrier copies.
   * With product support, the reader must return zero outside S and targets
   * are generated in its one-generator closure I. The adjoint returns the
   * full Hamiltonian image, including rows outside S; a selected-structure
   * caller projects those out in its exterior-transform adjoint.
   * A nonzero exchange_parity (+1 or -1) instead contracts one representative
   * of each alpha/beta exchange orbit with its exact multiplicity. This is
   * valid ONLY when both the input and the final left projection have that
   * transpose parity in identical spin-string spaces. The intermediate
   * output is then a projected-action representative, NOT the full H image.
   */
  void apply_tiled(
      int block_width,
      const std::function<void(const std::vector<int>&, const std::vector<int>&, Eigen::Ref<Eigen::MatrixXd>)>& read_columns,
      const std::function<void(const std::vector<int>&, const std::vector<int>&, const DirectCiColumnAccumulator&)>& add_columns,
      const ProductSupport* support = nullptr,
      int exchange_parity = 0) const;

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
   * A nonzero exchange_parity certifies X^T = parity * X for every input
   * block in identical alpha/beta string lists. Zero makes no such assumption.
   */
  DirectCiIntegralAdjoint integral_adjoint(
      const Eigen::Ref<const Eigen::MatrixXd>& left,
      const Eigen::Ref<const Eigen::MatrixXd>& right,
      int exchange_parity = 0) const;

  /** @brief Returns H * response and the integral derivative of response.dot(H * state).
   * Reuses the response scatter; two shared tiles serve both contractions.
   * Inputs follow apply() packing and may be non-owning strided blocks.
   * With certified exchange parity, sigma is reconstructed as a full image,
   * not the unprojected representative used by apply_tiled().
   */
  DirectCiResponse apply_with_adjoint(
      const Eigen::Ref<const Eigen::MatrixXd>& response,
      const Eigen::Ref<const Eigen::MatrixXd>& state,
      int exchange_parity = 0) const;

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

  /** @brief Shared channel sweep for sigma, integral adjoints, or both.
   * Non-null sigma receives H * left. A sigma-only request neither constructs
   * a kernel gradient nor scatters/copies the unused right-hand channels.
   */
  DirectCiIntegralAdjoint contract_channels(
      const Eigen::Ref<const Eigen::MatrixXd>& left,
      const Eigen::Ref<const Eigen::MatrixXd>& right,
      Eigen::MatrixXd* sigma, int exchange_parity, bool collect_integrals) const;

  /** @brief Cooperative sigma contraction inside the current OpenMP team.
   * Each output beta column has one owner. All workers must call this method.
   */
  void apply_sigma_tile(
      int alpha_begin, int alpha_count, int block_begin, int block_count,
      int beta_begin, int beta_count, int exchange_parity,
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
