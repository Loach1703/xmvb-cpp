#pragma once

#include <cstddef>
#include <functional>

#include <Eigen/Core>

#include "mcscf/ci/sigma.hpp"

namespace xmvb::mcscf {

/** @brief Exact orthogonal spin factors supplied a bounded source panel at a time.
 * load(first,count,A,B) fills the borrowed full-row views so that the represented
 * amplitude is sum_j A_j B_j^T. No complete product or all-source array is needed.
 */
struct SpinFactorSource {
  int n_alpha = 0, n_beta = 0, count = 0;
  std::function<void(int, int, Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXd>)> load;
};

/** @brief Cross-spin Hamiltonian factors after their same-spin use is complete.
 * alpha contains E_(alpha,P) A, beta contains sum_Q K_PQ E_(beta,Q) B.
 * Column P*width+j belongs to one channel/source. The consumer may overwrite
 * the views, but must finish before returning and may not retain them.
 */
struct FactorCrossBlock {
  int source, width, channel, channels;
  Eigen::Ref<Eigen::MatrixXd> alpha, beta;
};

/** @brief Base and same-spin images; all channel consumers have finished.
 * same_alpha/beta contain D_alpha A and D_beta B, or zero columns for an
 * overlap-only workspace. These borrowed views may be consumed in place.
 */
struct FactorImageBlock {
  int source, width;
  Eigen::Ref<Eigen::MatrixXd> alpha, beta, same_alpha, same_beta;
};

/** @brief Work and explicit numerical-slab bound, not total process memory. */
struct FactorActionWork {
  int source_tile = 0, channel_tile = 0;
  std::size_t panels = 0, sources = 0, channel_blocks = 0, workspace_bound = 0;
  double scatter_seconds = 0, mix_seconds = 0, gather_seconds = 0;
};

/** @brief Shared Slater--Condon source/channel sweep, reusable across RHS.
 * Owns one set of numerical slabs, shared by the OpenMP team. The workspace
 * budget controls only residency, never a physical-space approximation.
 */
class FactorActionWorkspace {
public:
  FactorActionWorkspace(const SpinFactorSource& source, int channels,
      std::size_t workspace_bytes, bool overlap_only = false);

  FactorActionWork apply(const DirectCiSigmaAction& sigma, const SpinFactorSource& source,
      const std::function<void(FactorCrossBlock)>& cross,
      const std::function<void(FactorImageBlock)>& finish);

private:
  int n_alpha_, n_beta_, channels_, source_tile_, channel_tile_;
  bool overlap_only_;
  Eigen::MatrixXd alpha_, beta_, ha_, hb_;
  Eigen::MatrixXd ap_, bp_, va_, vb_, aq_, bq_;
};

/** @brief Exact derivative of <left, H right> with respect to active integrals.
 * Sources describe the amplitudes of this bilinear form. Any parity quotient
 * must be justified by the caller's left projection; it is never inferred here.
 * Rank and channel panels are streamed; no alpha-by-beta array is allocated.
 * The budget covers numerical workspace, excluding immutable input metadata
 * and the required returned orbital-pair gradient. Left/right may differ.
 */
DirectCiIntegralAdjoint factor_integral_adjoint(const DirectCiSigmaAction& sigma,
    const SpinFactorSource& left, const SpinFactorSource& right,
    std::size_t workspace_bytes = 64ULL * 1024ULL * 1024ULL);

} // namespace xmvb::mcscf
