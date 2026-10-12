#pragma once

#include <vector>

#include <Eigen/Core>

#include "vbscf/integrals/ao/contracts/input.hpp"

namespace xmvb::vb {

/** @brief Shared inputs and reused worker outputs for two symmetric channels. */
struct AoH1eSymmetricWorkspace {
  using Channels = Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>;
  Channels sources;
  std::vector<Channels> partial;
};

/** @brief Reused interleaved batch inputs and AO-sized worker reductions. */
struct AoH1eBatchWorkspace {
  using Channels = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  Channels sources;
  std::vector<Channels> partial;
  std::vector<Eigen::VectorXd> coulomb;
};

/**
 * @brief Applies one physical symmetric response per source in a shared sweep.
 *
 * Columns are flattened symmetric AO matrices. Unlike the forward/adjoint
 * interface, each direction has only one channel and one complete output.
 * Sources and integral storage remain shared; worker reductions are reused.
 */
void apply_symmetric_ao_h1e_batch(
    const Eigen::Ref<const Eigen::MatrixXd>& sources,
    const AoIntegralInput& ao, int n_threads,
    AoH1eBatchWorkspace* workspace, Eigen::MatrixXd* response);

/** @brief Consume already interleaved workspace sources without copying them.
 * Each worker retains only packed symmetric AO output channels. The returned
 * columns are complete physical AO matrices, including doubled raw diagonals.
 */
void apply_symmetric_ao_h1e_channels(const AoIntegralInput& ao, int n_threads,
    AoH1eBatchWorkspace* workspace, Eigen::MatrixXd* response);

/**
 * @brief Applies two symmetric AO sources with a common two-lane stencil.
 *
 * Both inputs must be symmetric. The first output is the complete physical
 * AO forward matrix; the second is half the physical adjoint matrix, so its
 * sum with its transpose equals that of the raw unique-ERI transpose.
 * Input channels are shared; workers retain only AO-sized output buffers.
 */
void apply_symmetric_ao_h1e_fused(
    const double* source, const double* adjoint, const AoIntegralInput& ao,
    int n_threads, AoH1eSymmetricWorkspace* workspace,
    std::vector<double>* forward, std::vector<double>* transpose);

/** @brief Applies the exact AO-H1E operator directly from unique AO ERIs. */
std::vector<double> apply_ao_h1e(
    const double* source,
    const AoIntegralInput& ao,
    int n_threads);

/** @brief Applies the transpose exact AO-H1E operator directly from AO ERIs. */
std::vector<double> apply_ao_h1e_transpose(
    const double* adjoint,
    const AoIntegralInput& ao,
    int n_threads);

/**
 * @brief Applies symmetric sources to several directions in one ERI sweep.
 *
 * Columns are flattened symmetric AO matrices. Output conventions match
 * apply_symmetric_ao_h1e_fused: complete forward and half physical adjoint.
 * The integral stream and packed source channels are shared by all workers.
 */
void apply_symmetric_ao_h1e_fused_batch(
    const Eigen::Ref<const Eigen::MatrixXd>& sources,
    const Eigen::Ref<const Eigen::MatrixXd>& adjoints,
    const AoIntegralInput& ao,
    int n_threads,
    AoH1eBatchWorkspace* workspace,
    Eigen::MatrixXd* forward,
    Eigen::MatrixXd* transpose);

}  // namespace xmvb::vb
