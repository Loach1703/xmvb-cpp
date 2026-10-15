#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace xmvb::vb {

class StructureAction;
class ExteriorOrbitalTransform;

/** Complete-sector inverse used as an SPD selected-metric preconditioner.
 *
 * The first columns of the completed basis are the original selected VB
 * structures, with their original coefficients and ordering. The remaining
 * columns span only the missing spin-pure directions. No integral action or
 * dense structure overlap is built.
 */
class RestrictedMetricInverse {
public:
  class Prepared {
  public:
    /** Applies the principal full-inverse block N_PP to selected RHS columns.
     * The originating RestrictedMetricInverse must outlive this chart.
     * N_PP is SPD but is generally not the inverse of the selected Gram.
     */
    Eigen::MatrixXd precondition(const Eigen::Ref<const Eigen::MatrixXd>& rhs) const;
    int n_complement() const noexcept;
    std::size_t stored_bytes() const noexcept;

  private:
    friend class RestrictedMetricInverse;
    explicit Prepared(const RestrictedMetricInverse* owner);
    const RestrictedMetricInverse* owner_;
    std::shared_ptr<const ExteriorOrbitalTransform> alpha_;
    std::shared_ptr<const ExteriorOrbitalTransform> beta_;
  };

  /** Returns null when the carrier, spin sector, complement, or budget is unsuitable. */
  static std::unique_ptr<RestrictedMetricInverse> make(
      const StructureAction& action,
      const std::vector<std::vector<int>>& alpha_strings,
      const std::vector<std::vector<int>>& beta_strings,
      std::size_t byte_limit);

  /** Prepares one fixed-chart principal inverse action for SPD PCG. */
  std::shared_ptr<const Prepared> prepare(
      const Eigen::Ref<const Eigen::MatrixXd>& upper) const;

  int n_selected() const noexcept { return n_selected_; }
  int n_complement() const noexcept { return n_full_ - n_selected_; }
  std::size_t stored_bytes() const noexcept { return stored_bytes_; }

private:
  struct Block {
    std::vector<int> indices;
    Eigen::MatrixXd lower;
  };

  RestrictedMetricInverse(const StructureAction& action,
                          const std::vector<std::vector<int>>& alpha_strings,
                          const std::vector<std::vector<int>>& beta_strings,
                          std::size_t byte_limit);

  Eigen::MatrixXd apply_full(
      const Eigen::Ref<const Eigen::MatrixXd>& rhs,
      const ExteriorOrbitalTransform& alpha,
      const ExteriorOrbitalTransform& beta) const;
  Eigen::MatrixXd solve_spin(
      const Eigen::Ref<const Eigen::MatrixXd>& rhs) const;

  std::vector<std::vector<int>> alpha_strings_;
  std::vector<std::vector<int>> beta_strings_;
  Eigen::SparseMatrix<double, Eigen::ColMajor, int> expansion_;
  std::vector<Block> blocks_;
  int n_selected_ = 0;
  int n_full_ = 0;
  int n_alpha_ = 0;
  int n_beta_ = 0;
  int n_orbitals_ = 0;
  int max_block_size_ = 0;
  int chunk_width_ = 1;
  std::size_t stored_bytes_ = 0;
  std::size_t byte_limit_ = 0;
};

}  // namespace xmvb::vb
