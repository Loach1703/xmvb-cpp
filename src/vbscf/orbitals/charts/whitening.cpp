#include "vbscf/orbitals/charts/whitening.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <Eigen/LU>
#include <lapacke.h>

#include "vbscf/diagnostics/performance/phases.hpp"

namespace xmvb::vb {

WhitenedTangent whiten_orbital_tangent(
    const Eigen::Ref<const Eigen::MatrixXd>& quotient,
    const Eigen::Ref<const Eigen::VectorXd>& coefficients,
    const Eigen::Ref<const Eigen::MatrixXd>& overlap,
    bool sphere,
    const Eigen::MatrixXd* quotient_metric) {
  diagnostics::ScopedPhase phase(diagnostics::Phase::ChartWhitening);
  if (quotient.rows() != coefficients.size() ||
      overlap.rows() != coefficients.size() || overlap.cols() != coefficients.size())
    throw std::invalid_argument("normalized-orbital quotient metric dimensions are inconsistent");
  if (!quotient.cols()) return {quotient, Eigen::MatrixXd(0, 0)};

  const Eigen::VectorXd overlap_raw = overlap * coefficients;
  const double squared_norm = coefficients.dot(overlap_raw);
  if (!(squared_norm > 0.0) || !std::isfinite(squared_norm))
    throw std::runtime_error("normalized-orbital quotient metric has nonpositive orbital norm");
  const double norm = std::sqrt(squared_norm);
  const Eigen::VectorXd normalized = coefficients / norm;
  const Eigen::VectorXd radial = quotient.transpose() * (overlap_raw / norm);
  if (quotient_metric && (quotient_metric->rows() != quotient.cols() ||
                         quotient_metric->cols() != quotient.cols()))
    throw std::invalid_argument("shared quotient metric dimensions are inconsistent");
  // JQ=(Q-c r^T)/||b||_S and c^T S c=1, hence
  // G=((Q^T S Q)-r r^T)/||b||_S^2. Q^T S Q can be shared by full-AO
  // orbitals with the same accepted-point gauge; each radial update remains exact.
  Eigen::MatrixXd metric = quotient_metric ? *quotient_metric
      : (quotient.transpose() * overlap * quotient).eval();
  metric.noalias() -= radial * radial.transpose();
  metric /= squared_norm;
  metric = (0.5 * (metric + metric.transpose())).eval();

  // Retain the unique symmetric inverse square root used by the chart.
  // LAPACK changes the implementation, not the coordinate convention.
  Eigen::MatrixXd eigenvectors = metric;
  Eigen::VectorXd eigenvalues(metric.rows());
  const lapack_int dimension = static_cast<lapack_int>(metric.rows());
  const lapack_int info = LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'L', dimension,
      eigenvectors.data(), dimension, eigenvalues.data());
  if (info != 0 || !eigenvalues.allFinite() || !eigenvectors.allFinite())
    throw std::runtime_error("failed to factor normalized-orbital quotient metric");
  const double spectral_scale = eigenvalues.maxCoeff();
  const double floor = 64.0 * std::numeric_limits<double>::epsilon() *
      static_cast<double>(metric.rows()) * spectral_scale;
  if (!(spectral_scale > 0.0) || eigenvalues.minCoeff() <= floor)
    throw std::runtime_error("normalized-orbital quotient metric is rank deficient");

  const Eigen::MatrixXd inverse_root = eigenvectors *
      eigenvalues.cwiseSqrt().cwiseInverse().asDiagonal() * eigenvectors.transpose();
  WhitenedTangent result;
  if (sphere) {
    Eigen::MatrixXd physical = quotient;
    physical.noalias() -= normalized * radial.transpose();
    result.basis.noalias() = physical * inverse_root;
  } else {
    result.basis.noalias() = quotient * inverse_root;
  }
  diagnostics::ScopedPhase gram_phase(diagnostics::Phase::ChartGram);
  result.inverse_gram = metric;
  if (sphere) {
    // Q^T Q = I. With c normalized and r=Q^T S c, the raw tangent is
    // (Q-c r^T) G^{-1/2}. Its unwhitened Gram is I+Z D Z^T, where
    // Z=[Q^T c,r], D=[[0,-1],[-1,c^T c]]. Apply the 2x2 Woodbury identity;
    // this is exact even when Q^T c has a finite numerical residual.
    Eigen::MatrixXd update(metric.rows(), 2);
    update.col(0).noalias() = quotient.transpose() * normalized;
    update.col(1) = radial;
    Eigen::Matrix2d small = update.transpose() * update;
    small(0, 0) -= normalized.squaredNorm();
    small(0, 1) -= 1.0;
    small(1, 0) -= 1.0;
    const Eigen::FullPivLU<Eigen::Matrix2d> factor(small);
    if (!factor.isInvertible())
      throw std::runtime_error("singular normalized-orbital raw tangent metric");
    const Eigen::MatrixXd weighted = eigenvectors *
        eigenvalues.cwiseSqrt().asDiagonal() * (eigenvectors.transpose() * update);
    result.inverse_gram.noalias() -= weighted * factor.solve(weighted.transpose());
    result.inverse_gram = (0.5 * (result.inverse_gram + result.inverse_gram.transpose())).eval();
  }
  if (!result.basis.allFinite() || !result.inverse_gram.allFinite())
    throw std::runtime_error("non-finite normalized-orbital tangent factor");
  return result;
}

}  // namespace xmvb::vb
