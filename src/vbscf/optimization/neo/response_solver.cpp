#include "vbscf/optimization/neo/response_solver.hpp"
#include "vbscf/diagnostics/performance/phases.hpp"
#include "vbscf/diagnostics/performance/neo_trace.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <lapacke.h>

#include "vbscf/optimization/trust_region/spectral.hpp"

namespace xmvb::vb {
namespace {

bool orthonormalize(
    Eigen::VectorXd* direction,
    const Eigen::Ref<const Eigen::MatrixXd>& basis,
    double* retained_fraction = nullptr) {
  if (retained_fraction) *retained_fraction = 0.0;
  const double initial_norm = direction->stableNorm();
  if (!(initial_norm > 0.0) || !std::isfinite(initial_norm)) return false;
  for (int pass = 0; pass < 2; ++pass) {
    direction->noalias() -= basis * (basis.transpose() * *direction);
  }
  const double norm = direction->stableNorm();
  if (retained_fraction) *retained_fraction = norm / initial_norm;
  const double threshold = 64.0 * std::numeric_limits<double>::epsilon() *
      initial_norm;
  if (!(norm > threshold) || !std::isfinite(norm)) return false;
  *direction /= norm;
  return true;
}

void append_column(Eigen::MatrixXd* matrix, const Eigen::VectorXd& column) {
  if (matrix->cols() == 0) {
    matrix->resize(column.size(), 1);
  } else {
    matrix->conservativeResize(Eigen::NoChange, matrix->cols() + 1);
  }
  matrix->col(matrix->cols() - 1) = column;
}

std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

Eigen::VectorXd generic_probe(Eigen::Index size) {
  Eigen::VectorXd probe(size);
  for (Eigen::Index j = 0; j < size; ++j) {
    const std::uint64_t bits = splitmix64(static_cast<std::uint64_t>(j));
    const double unit = static_cast<double>(bits >> 11U) * 0x1.0p-53;
    probe[j] = 2.0 * unit - 1.0;
  }
  return probe;
}

Eigen::MatrixXd symmetric_part(
    const Eigen::MatrixXd& matrix,
    const char* message,
    double operator_relative_accuracy = 0.0) {
  const double factor = std::max(
      2048.0 * std::numeric_limits<double>::epsilon() *
          static_cast<double>(std::max<Eigen::Index>(1, matrix.rows())),
      operator_relative_accuracy);
  const double defect = (matrix - matrix.transpose()).stableNorm();
  const double scale = std::max(1.0, matrix.stableNorm());
  if (defect > factor * scale) {
    std::ostringstream detail;
    detail.precision(17);
    detail << message << ": relative defect=" << defect / scale;
    throw std::runtime_error(detail.str());
  }
  return 0.5 * (matrix + matrix.transpose());
}

struct StructureInverse {
  Eigen::MatrixXd value;
  Eigen::MatrixXd null_vectors;
};

StructureInverse symmetric_pseudoinverse(const Eigen::MatrixXd& matrix) {
  if (matrix.rows() == 0) {
    return {matrix, Eigen::MatrixXd(matrix.rows(), 0)};
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(matrix);
  if (spectrum.info() != Eigen::Success) {
    throw std::runtime_error("NEO structure eigensystem failed");
  }
  const double cutoff = 256.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, spectrum.eigenvalues().cwiseAbs().maxCoeff()) *
      static_cast<double>(matrix.rows());
  Eigen::VectorXd inverse = spectrum.eigenvalues();
  Eigen::Index nullity = 0;
  for (Eigen::Index j = 0; j < inverse.size(); ++j) {
    if (std::abs(inverse[j]) > cutoff) {
      inverse[j] = 1.0 / inverse[j];
    } else {
      inverse[j] = 0.0;
      ++nullity;
    }
  }
  StructureInverse result;
  result.value = spectrum.eigenvectors() * inverse.asDiagonal() *
      spectrum.eigenvectors().transpose();
  result.null_vectors.resize(matrix.rows(), nullity);
  Eigen::Index column = 0;
  for (Eigen::Index j = 0; j < inverse.size(); ++j) {
    if (inverse[j] == 0.0) {
      result.null_vectors.col(column++) = spectrum.eigenvectors().col(j);
    }
  }
  return result;
}

Eigen::VectorXd preconditioned_orbital(
    const ResponseNeoProblem& problem,
    const Eigen::VectorXd& residual,
    double relative_residual_target,
    double shift = 0.0) {
  return problem.has_orbital_preconditioner()
      ? problem.apply_orbital_preconditioner(
            residual, shift, relative_residual_target)
      : residual;
}

double relative_preconditioner_budget(
    double absolute_target,
    double covector_norm) {
  return covector_norm == 0.0
      ? 1.0
      : std::min(1.0, absolute_target / covector_norm);
}

}  // namespace

ResponseNeoProblem::ResponseNeoProblem(
    Eigen::VectorXd orbital_gradient,
    Eigen::Index structure_size,
    ResponseNeoBlockAction apply_orbital_coupling,
    ResponseNeoBlockAction apply_structure_coupling,
    NeoAction apply_orbital_metric,
    NeoDualNorm orbital_dual_norm,
    ResponseNeoPreconditioner apply_orbital_preconditioner,
    Eigen::VectorXd initial_orbital_guess,
    double operator_relative_accuracy,
    ResponseNeoMatrixAction apply_structure_coupling_block,
    ResponseNeoStructureResidualSolver solve_structure_residual)
    : orbital_gradient_(std::move(orbital_gradient)),
      structure_size_(structure_size),
      apply_orbital_coupling_(std::move(apply_orbital_coupling)),
      apply_structure_coupling_(std::move(apply_structure_coupling)),
      apply_structure_coupling_block_(
          std::move(apply_structure_coupling_block)),
      apply_orbital_metric_(std::move(apply_orbital_metric)),
      orbital_dual_norm_(std::move(orbital_dual_norm)),
      apply_orbital_preconditioner_(std::move(apply_orbital_preconditioner)),
      solve_structure_residual_(std::move(solve_structure_residual)),
      initial_orbital_guess_(std::move(initial_orbital_guess)),
      operator_relative_accuracy_(operator_relative_accuracy) {
  if (orbital_gradient_.size() == 0 || !orbital_gradient_.allFinite() ||
      structure_size_ < 0) {
    throw std::invalid_argument("response NEO requires finite dimensions");
  }
  if (!apply_orbital_coupling_ || !apply_structure_coupling_ ||
      !apply_orbital_metric_ || !orbital_dual_norm_) {
    throw std::invalid_argument(
        "response NEO requires coupling, metric, and dual-norm actions");
  }
  if (initial_orbital_guess_.size() != 0 &&
      (initial_orbital_guess_.size() != orbital_size() ||
       !initial_orbital_guess_.allFinite())) {
    throw std::invalid_argument(
        "response NEO initial orbital guess is invalid");
  }
  if (!std::isfinite(operator_relative_accuracy_) ||
      operator_relative_accuracy_ < 0.0) {
    throw std::invalid_argument(
        "response NEO operator accuracy must be finite and nonnegative");
  }
}

ResponseNeoDirection ResponseNeoProblem::apply_block_checked(
    const ResponseNeoBlockAction& action,
    const Eigen::VectorXd& vector,
    Eigen::Index expected_size) const {
  if (vector.size() != expected_size || !vector.allFinite()) {
    throw std::invalid_argument("invalid vector passed to NEO block action");
  }
  ResponseNeoDirection result = action(vector);
  if (result.orbital.size() != orbital_size() ||
      result.structure.size() != structure_size() ||
      !result.orbital.allFinite() || !result.structure.allFinite()) {
    throw std::runtime_error("coupled NEO action returned an invalid vector");
  }
  return result;
}

ResponseNeoDirection ResponseNeoProblem::apply_orbital_coupling(
    const Eigen::VectorXd& direction) const {
  return apply_block_checked(
      apply_orbital_coupling_, direction, orbital_size());
}

ResponseNeoDirection ResponseNeoProblem::apply_structure_coupling(
    const Eigen::VectorXd& direction) const {
  return apply_block_checked(
      apply_structure_coupling_, direction, structure_size());
}

ResponseNeoDirectionBlock ResponseNeoProblem::apply_structure_coupling_block(
    const Eigen::MatrixXd& directions) const {
  if (directions.rows() != structure_size() || directions.cols() == 0 ||
      !directions.allFinite()) {
    throw std::invalid_argument(
        "invalid matrix passed to NEO structure block action");
  }
  ResponseNeoDirectionBlock result;
  if (apply_structure_coupling_block_) {
    result = apply_structure_coupling_block_(directions);
  } else {
    result.orbital.resize(orbital_size(), directions.cols());
    result.structure.resize(structure_size(), directions.cols());
    for (Eigen::Index column = 0; column < directions.cols(); ++column) {
      const ResponseNeoDirection image =
          apply_structure_coupling(directions.col(column));
      result.orbital.col(column) = image.orbital;
      result.structure.col(column) = image.structure;
    }
  }
  if (result.orbital.rows() != orbital_size() ||
      result.structure.rows() != structure_size() ||
      result.orbital.cols() != directions.cols() ||
      result.structure.cols() != directions.cols() ||
      !result.orbital.allFinite() || !result.structure.allFinite()) {
    throw std::runtime_error(
        "coupled NEO block action returned an invalid matrix");
  }
  return result;
}

Eigen::VectorXd ResponseNeoProblem::apply_orbital_checked(
    const NeoAction& action,
    const Eigen::VectorXd& vector,
    const char* message) const {
  if (vector.size() != orbital_size() || !vector.allFinite()) {
    throw std::invalid_argument("invalid orbital vector passed to response NEO");
  }
  Eigen::VectorXd result = action(vector);
  if (result.size() != orbital_size() || !result.allFinite()) {
    throw std::runtime_error(message);
  }
  return result;
}

Eigen::VectorXd ResponseNeoProblem::apply_orbital_metric(
    const Eigen::VectorXd& direction) const {
  diagnostics::ScopedPhase phase(diagnostics::Phase::NeoMetric);
  return apply_orbital_checked(
      apply_orbital_metric_, direction,
      "response NEO orbital metric returned an invalid vector");
}

double ResponseNeoProblem::orbital_dual_norm(
    const Eigen::VectorXd& covector) const {
  diagnostics::ScopedPhase phase(diagnostics::Phase::NeoMetric);
  if (covector.size() != orbital_size() || !covector.allFinite()) {
    throw std::invalid_argument(
        "invalid orbital covector passed to response NEO dual norm");
  }
  const double norm = orbital_dual_norm_(covector);
  if (!std::isfinite(norm) || norm < 0.0) {
    throw std::runtime_error(
        "response NEO orbital dual norm returned an invalid value");
  }
  return norm;
}

Eigen::VectorXd ResponseNeoProblem::apply_orbital_preconditioner(
    const Eigen::VectorXd& covector,
    double shift,
    double relative_residual_target) const {
  diagnostics::ScopedPhase phase(diagnostics::Phase::NeoPreconditioner);
  if (!apply_orbital_preconditioner_) {
    throw std::logic_error("response NEO has no orbital preconditioner");
  }
  if (!(shift >= 0.0) || !std::isfinite(shift)) {
    throw std::invalid_argument(
        "response NEO preconditioner shift must be finite and nonnegative");
  }
  if (!(relative_residual_target >= 0.0) ||
      !std::isfinite(relative_residual_target)) {
    throw std::invalid_argument(
        "response NEO preconditioner fraction must be finite and nonnegative");
  }
  if (covector.size() != orbital_size() || !covector.allFinite()) {
    throw std::invalid_argument(
        "invalid vector passed to response NEO preconditioner");
  }
  Eigen::VectorXd result = apply_orbital_preconditioner_(
      covector, shift, relative_residual_target);
  if (result.size() != orbital_size() || !result.allFinite()) {
    throw std::runtime_error(
        "response NEO orbital preconditioner returned an invalid vector");
  }
  return result;
}

Eigen::VectorXd ResponseNeoProblem::solve_structure_residual(
    const Eigen::VectorXd& residual) const {
  if (!solve_structure_residual_) {
    throw std::logic_error("response NEO has no structure residual solver");
  }
  if (residual.size() != structure_size() || !residual.allFinite()) {
    throw std::invalid_argument(
        "invalid vector passed to response NEO structure residual solver");
  }
  Eigen::VectorXd correction = solve_structure_residual_(residual);
  if (correction.size() != structure_size() || !correction.allFinite()) {
    throw std::runtime_error(
        "response NEO structure residual solver returned an invalid vector");
  }
  return correction;
}

bool ResponseNeoWorkspace::append_orbital(
    Eigen::VectorXd direction,
    int* actions,
    int* orbital_actions) {
  diagnostics::ScopedPhase phase(diagnostics::Phase::NeoBasis);
  double retained_fraction = 0.0;
  const bool independent = orthonormalize(
      &direction, orbital_basis_,
      diagnostics::neo_trace_enabled() ? &retained_fraction : nullptr);
  diagnostics::neo_trace("basis", "source=%s,independent=%d,retained=%.17g,size=%ld",
      diagnostics::neo_trace_context.source, static_cast<int>(independent),
      retained_fraction, static_cast<long>(orbital_basis_.cols()));
  if (!independent) return false;
  const Eigen::VectorXd metric = problem_.apply_orbital_metric(direction);
  const ResponseNeoDirection image =
      problem_.apply_orbital_coupling(direction);
  append_column(&orbital_basis_, direction);
  append_column(&orbital_metric_images_, metric);
  append_column(&orbital_orbital_images_, image.orbital);
  append_column(&orbital_structure_images_, image.structure);
  const Eigen::Index n_orbital = orbital_basis_.cols();
  const Eigen::Index n_structure = structure_basis_.cols();
  projected_orbital_.conservativeResize(n_orbital, n_orbital);
  projected_orbital_.col(n_orbital - 1).noalias() =
      orbital_basis_.transpose() * image.orbital;
  projected_orbital_.row(n_orbital - 1).noalias() =
      direction.transpose() * orbital_orbital_images_;
  projected_metric_.conservativeResize(n_orbital, n_orbital);
  projected_metric_.col(n_orbital - 1).noalias() =
      orbital_basis_.transpose() * metric;
  projected_metric_.row(n_orbital - 1).noalias() =
      direction.transpose() * orbital_metric_images_;
  projected_coupling_.conservativeResize(n_structure, n_orbital);
  if (n_structure != 0) {
    projected_coupling_.col(n_orbital - 1).noalias() =
        structure_basis_.transpose() * image.structure;
  }
  projected_coupling_adjoint_.conservativeResize(n_orbital, n_structure);
  if (n_structure != 0) {
    projected_coupling_adjoint_.row(n_orbital - 1).noalias() =
        direction.transpose() * structure_orbital_images_;
  }
  ++*actions;
  ++*orbital_actions;
  return true;
}

bool ResponseNeoWorkspace::append_structure(
    Eigen::VectorXd direction,
    int* actions,
    int* structure_actions) {
  Eigen::MatrixXd directions(direction.size(), 1);
  directions.col(0) = std::move(direction);
  return append_structure_block(
      std::move(directions), actions, structure_actions);
}

void ResponseNeoWorkspace::store_structure(
    const Eigen::VectorXd& direction,
    const ResponseNeoDirection& image) {
  append_column(&structure_basis_, direction);
  append_column(&structure_orbital_images_, image.orbital);
  append_column(&structure_structure_images_, image.structure);
  const Eigen::Index n_orbital = orbital_basis_.cols();
  const Eigen::Index n_structure = structure_basis_.cols();
  projected_coupling_.conservativeResize(n_structure, n_orbital);
  if (n_orbital != 0) {
    projected_coupling_.row(n_structure - 1).noalias() =
        direction.transpose() * orbital_structure_images_;
  }
  projected_coupling_adjoint_.conservativeResize(n_orbital, n_structure);
  if (n_orbital != 0) {
    projected_coupling_adjoint_.col(n_structure - 1).noalias() =
        orbital_basis_.transpose() * image.orbital;
  }
  projected_structure_.conservativeResize(n_structure, n_structure);
  projected_structure_.col(n_structure - 1).noalias() =
      structure_basis_.transpose() * image.structure;
  projected_structure_.row(n_structure - 1).noalias() =
      direction.transpose() * structure_structure_images_;
}

bool ResponseNeoWorkspace::retain_structure_image(
    Eigen::VectorXd direction,
    ResponseNeoDirection image) {
  const double initial_norm = direction.stableNorm();
  if (!(initial_norm > 0.0) || !std::isfinite(initial_norm)) return false;
  // Apply the same linear changes to the vector and its coupled image.
  // H(q - Z a) = Hq - (HZ)a, so no new action is needed.
  for (int pass = 0; pass < 2; ++pass) {
    const Eigen::VectorXd coefficients =
        structure_basis_.transpose() * direction;
    direction.noalias() -= structure_basis_ * coefficients;
    image.orbital.noalias() -= structure_orbital_images_ * coefficients;
    image.structure.noalias() -= structure_structure_images_ * coefficients;
  }
  const double norm = direction.stableNorm();
  const double threshold =
      64.0 * std::numeric_limits<double>::epsilon() * initial_norm;
  if (!(norm > threshold) || !std::isfinite(norm)) return false;
  direction /= norm;
  image.orbital /= norm;
  image.structure /= norm;
  store_structure(direction, image);
  return true;
}

bool ResponseNeoWorkspace::append_structure_block(
    Eigen::MatrixXd directions,
    int* actions,
    int* structure_actions) {
  Eigen::MatrixXd accepted(problem_.structure_size(), 0);
  for (Eigen::Index column = 0; column < directions.cols(); ++column) {
    Eigen::VectorXd direction = directions.col(column);
    if (!orthonormalize(&direction, structure_basis_) ||
        !orthonormalize(&direction, accepted)) {
      continue;
    }
    append_column(&accepted, direction);
  }
  if (accepted.cols() == 0) return false;
  const ResponseNeoDirectionBlock images =
      problem_.apply_structure_coupling_block(accepted);
  for (Eigen::Index column = 0; column < accepted.cols(); ++column) {
    store_structure(
        accepted.col(column),
        ResponseNeoDirection{
            images.orbital.col(column), images.structure.col(column)});
  }
  *actions += static_cast<int>(accepted.cols());
  *structure_actions += static_cast<int>(accepted.cols());
  return true;
}

ResponseNeoResult ResponseNeoWorkspace::solve(
    const NeoOptions& options, bool allow_early_trial,
    const std::function<void(const ResponseNeoResult&)>& observe_candidate,
    const ResponseNeoModelObserver& observe_model) {
  diagnostics::ScopedPhase phase(diagnostics::Phase::NeoSolve);
  if (!(options.trust_radius > 0.0) ||
      !std::isfinite(options.trust_radius) ||
      !(options.relative_residual_tolerance > 0.0) ||
      !std::isfinite(options.relative_residual_tolerance) ||
      options.absolute_residual_tolerance < 0.0 ||
      !std::isfinite(options.absolute_residual_tolerance) ||
      options.maximum_subspace_dimension < 0) {
    throw std::invalid_argument("invalid response NEO solver options");
  }

  const int full_dimension = static_cast<int>(
      problem_.orbital_size() + problem_.structure_size());
  const int maximum_dimension = options.maximum_subspace_dimension == 0
      ? full_dimension
      : std::min(full_dimension, options.maximum_subspace_dimension);
  if (gradient_norm_ < 0.0) {
    gradient_norm_ = problem_.orbital_dual_norm(problem_.orbital_gradient());
  }
  const double gradient_norm = gradient_norm_;
  const double kkt_residual_target =
      options.absolute_residual_tolerance +
      options.relative_residual_tolerance * gradient_norm;
  const double numerical_curvature_tolerance =
      std::sqrt(std::numeric_limits<double>::epsilon()) *
      static_cast<double>(std::max(1, full_dimension));
  const double curvature_relative_tolerance = std::min(
      options.relative_residual_tolerance,
      numerical_curvature_tolerance);
  ResponseNeoResult result;
  result.step.orbital = Eigen::VectorXd::Zero(problem_.orbital_size());
  result.step.structure = Eigen::VectorXd::Zero(problem_.structure_size());
  result.hessian_step = result.step;
  result.orbital_metric_step = result.step.orbital;
  result.kkt_residual = result.step;
  result.minimum_curvature_orbital = result.step.orbital;

  if (orbital_basis_.cols() == 0) {
    diagnostics::neo_trace_context.source = "initial";
    Eigen::VectorXd first =
        problem_.initial_orbital_guess().size() == problem_.orbital_size()
            ? problem_.initial_orbital_guess()
            : preconditioned_orbital(
                  problem_, -problem_.orbital_gradient(),
                  relative_preconditioner_budget(
                      kkt_residual_target, gradient_norm));
    if (!append_orbital(
            std::move(first), &result.coupled_actions,
            &result.orbital_actions)) {
      append_orbital(
          generic_probe(problem_.orbital_size()), &result.coupled_actions,
          &result.orbital_actions);
    }
  }
  if (orbital_basis_.cols() == 0) {
    throw std::runtime_error("response NEO could not construct an orbital basis");
  }
  if (options.require_curvature_certificate &&
      orbital_basis_.cols() == 1 && maximum_dimension > 1 &&
      problem_.orbital_size() > 1) {
    append_orbital(
        generic_probe(problem_.orbital_size()),
        &result.coupled_actions,
        &result.orbital_actions);
  }
  while (true) {
    diagnostics::ScopedPhase projected_phase(diagnostics::Phase::NeoProjectedSolve);
    diagnostics::ScopedPhase factor_phase(diagnostics::Phase::NeoProjectionFactor);
    ++result.iterations;
    const Eigen::MatrixXd& qo = orbital_basis_;
    const Eigen::MatrixXd& mo = orbital_metric_images_;
    const Eigen::MatrixXd& ao = orbital_orbital_images_;
    const Eigen::MatrixXd& bo = orbital_structure_images_;
    const Eigen::MatrixXd& qs = structure_basis_;
    const Eigen::MatrixXd& bt = structure_orbital_images_;
    const Eigen::MatrixXd& cs = structure_structure_images_;

    if (model_orbital_size_ != qo.cols() ||
        model_structure_size_ != qs.cols()) {
      const Eigen::MatrixXd a = symmetric_part(
          projected_orbital_, "response NEO orbital block is not symmetric",
          problem_.operator_relative_accuracy());
      if (metric_size_ != qo.cols()) {
        diagnostics::record_matrix_shape(diagnostics::Phase::NeoMetricFactor,
                                        qo.cols(), qo.cols());
        diagnostics::ScopedPhase metric_phase(diagnostics::Phase::NeoMetricFactor);
        const Eigen::MatrixXd metric = symmetric_part(
            projected_metric_, "response NEO orbital metric is not symmetric");
        // The retained metric block is immutable. Extend R in M=R^T R and
        // its inverse instead of refactoring every leading principal block.
        const Eigen::Index previous = std::max<Eigen::Index>(0, metric_size_);
        metric_factor_.conservativeResize(metric.rows(), metric.cols());
        whitening_.conservativeResize(metric.rows(), metric.cols());
        metric_factor_.bottomRows(metric.rows() - previous).setZero();
        whitening_.bottomRows(metric.rows() - previous).setZero();
        for (Eigen::Index column = previous; column < metric.cols(); ++column) {
          const Eigen::VectorXd border = metric_factor_.topLeftCorner(column, column)
              .transpose().triangularView<Eigen::Lower>()
              .solve(metric.col(column).head(column));
          const double pivot = metric(column, column) - border.squaredNorm();
          if (!(pivot > 0.0) || !std::isfinite(pivot))
            throw std::runtime_error("response NEO orbital metric is not positive definite");
          const double diagonal = std::sqrt(pivot);
          metric_factor_.col(column).head(column) = border;
          metric_factor_(column, column) = diagonal;
          whitening_.col(column).head(column).noalias() =
              -whitening_.topLeftCorner(column, column) * border / diagonal;
          whitening_(column, column) = 1.0 / diagonal;
        }
        metric_size_ = qo.cols();
      }

      Eigen::MatrixXd b(qs.cols(), qo.cols());
      if (inverse_size_ != qs.cols()) {
        diagnostics::record_matrix_shape(diagnostics::Phase::NeoStructureFactor,
                                        qs.cols(), qs.cols());
        diagnostics::ScopedPhase structure_phase(diagnostics::Phase::NeoStructureFactor);
        const Eigen::MatrixXd c = symmetric_part(
            projected_structure_,
            "response NEO structure block is not symmetric",
            problem_.operator_relative_accuracy());
        StructureInverse inverse = symmetric_pseudoinverse(c);
        structure_inverse_ = std::move(inverse.value);
        structure_null_vectors_ = std::move(inverse.null_vectors);
        inverse_size_ = qs.cols();
      }
      if (qs.cols() == 0) {
        b.setZero();
      } else {
        b = projected_coupling_;
        const Eigen::MatrixXd& projected_adjoint =
            projected_coupling_adjoint_;
        const double adjoint_scale = std::max({
            1.0, b.stableNorm(), projected_adjoint.stableNorm()});
        const double adjoint_tolerance = std::max(
            256.0 * std::numeric_limits<double>::epsilon() *
                static_cast<double>(std::max(b.rows(), b.cols())),
            problem_.operator_relative_accuracy());
        if ((b.transpose() - projected_adjoint).stableNorm() >
            adjoint_tolerance * adjoint_scale) {
          throw std::runtime_error("response NEO coupling blocks are not adjoints");
        }
        b = 0.5 * (b + projected_adjoint.transpose());
        if (structure_null_vectors_.cols() != 0) {
          const Eigen::MatrixXd range_defect =
              structure_null_vectors_.transpose() * b;
          const double range_error = range_defect.stableNorm();
          const double range_tolerance = std::max(
              256.0 * std::numeric_limits<double>::epsilon() *
                  static_cast<double>(std::max(qs.cols(), b.cols())),
              problem_.operator_relative_accuracy()) *
              std::max(1.0, b.stableNorm());
          if (range_error > range_tolerance) {
            Eigen::Index unresolved = 0;
            range_defect.rowwise().squaredNorm().maxCoeff(&unresolved);
            if (structure_basis_.cols() == problem_.structure_size()) {
              throw std::runtime_error(
                  "response NEO coupling is outside the range of C");
            }
            const int basis_size = static_cast<int>(
                orbital_basis_.cols() + structure_basis_.cols());
            if (basis_size >= maximum_dimension) {
              result.stop_reason = NeoStopReason::SubspaceLimit;
              return result;
            }
            if (append_structure(
                    cs * structure_null_vectors_.col(unresolved),
                    &result.coupled_actions,
                    &result.structure_actions)) {
              continue;
            }
            throw std::runtime_error(
                "response NEO coupling is outside the range of C");
          }
        }
      }

      response_coefficients_.noalias() = -structure_inverse_ * b;
      const Eigen::MatrixXd relaxed = a + b.transpose() * response_coefficients_;
      if (observe_model) observe_model(qo, a, relaxed);
      spectrum_vectors_ =
          whitening_.transpose() * relaxed * whitening_;
      diagnostics::ScopedPhase spectrum_phase(diagnostics::Phase::NeoSpectrum);
      diagnostics::record_matrix_shape(diagnostics::Phase::NeoSpectrum,
                                      qo.cols(), qo.cols());
      const lapack_int dimension = static_cast<lapack_int>(qo.cols());
      spectrum_values_.resize(dimension);
      // The projected model needs its complete spectrum. Divide-and-conquer
      // computes that spectrum in the owned small matrix, avoiding Eigen's
      // QR eigensolve followed by a second eigenvector-matrix copy.
      const lapack_int info = LAPACKE_dsyevd(
          LAPACK_COL_MAJOR, 'V', 'L', dimension,
          spectrum_vectors_.data(), dimension, spectrum_values_.data());
      if (info != 0) {
        result.stop_reason = NeoStopReason::NumericalFailure;
        return result;
      }
      const Eigen::VectorXd gradient =
          qo.transpose() * problem_.orbital_gradient();
      const Eigen::VectorXd white_gradient = whitening_.transpose() * gradient;
      spectral_gradient_.noalias() = spectrum_vectors_.transpose() * white_gradient;
      model_orbital_size_ = qo.cols();
      model_structure_size_ = qs.cols();
    }
    const Eigen::MatrixXd& whitening = whitening_;
    const Eigen::MatrixXd& response_coefficients = response_coefficients_;
    const SpectralTrustRegionSolution projected = solve_spectral_trust_region(
        spectrum_values_, spectral_gradient_,
        options.trust_radius);
    const Eigen::VectorXd p_coefficients = whitening *
        spectrum_vectors_ * projected.step;
    const Eigen::VectorXd q_coefficients =
        response_coefficients * p_coefficients;

    factor_phase.stop();
    diagnostics::ScopedPhase residual_phase(diagnostics::Phase::NeoProjectionResidual);

    result.step.orbital.noalias() = qo * p_coefficients;
    result.step.structure.noalias() = qs * q_coefficients;
    const Eigen::VectorXd orbital_a = ao * p_coefficients;
    const Eigen::VectorXd orbital_bt = bt * q_coefficients;
    const Eigen::VectorXd structure_b = bo * p_coefficients;
    const Eigen::VectorXd structure_c = cs * q_coefficients;
    result.hessian_step.orbital = orbital_a + orbital_bt;
    result.hessian_step.structure = structure_b + structure_c;
    result.orbital_metric_step.noalias() = mo * p_coefficients;
    result.shift = projected.shift;
    result.boundary = projected.boundary;
    result.hard_case = projected.hard_case;
    result.step_norm = std::sqrt(std::max(
        0.0, result.step.orbital.dot(result.orbital_metric_step)));
    result.kkt_residual.orbital = problem_.orbital_gradient() +
        result.hessian_step.orbital +
        result.shift * result.orbital_metric_step;
    result.kkt_residual.structure = result.hessian_step.structure;
    // Inexact Newton convergence is controlled relative to the current
    // gradient, not to the much larger cancelling Hessian terms.  Scaling by
    // those terms would permit an O(||g||) residual in the physical dual
    // metric and destroy the local Newton rate in ill-conditioned coordinates.
    const double orbital_residual_norm =
        problem_.orbital_dual_norm(result.kkt_residual.orbital);
    const double structure_residual_norm =
        result.kkt_residual.structure.stableNorm();
    result.residual_norm = std::hypot(
        orbital_residual_norm, structure_residual_norm);
    result.residual_target = kkt_residual_target;
    if (diagnostics::neo_trace_enabled() &&
        ((result.orbital_actions & (result.orbital_actions - 1)) == 0 ||
         result.residual_norm <= kkt_residual_target)) {
      // Reuse retained action images. These coordinate norms diagnose model
      // cancellation; they never replace the physical KKT certificate.
      const Eigen::VectorXd fixed_residual = problem_.orbital_gradient() +
          orbital_a + result.shift * result.orbital_metric_step;
      const double fixed_norm = fixed_residual.stableNorm();
      const double response_norm = orbital_bt.stableNorm();
      const double product = fixed_norm * response_norm;
      diagnostics::neo_trace("residual_terms",
          "hvp=%d,no=%ld,ns=%ld,gradient=%.17g,fixed=%.17g,response=%.17g,cosine=%.17g,combined=%.17g,forcing=%.17g,structure_image=%.17g",
          result.orbital_actions, static_cast<long>(qo.cols()),
          static_cast<long>(qs.cols()), problem_.orbital_gradient().stableNorm(),
          fixed_norm, response_norm,
          product > 0.0 ? fixed_residual.dot(orbital_bt) / product : 0.0,
          result.kkt_residual.orbital.stableNorm(),
          structure_b.stableNorm(), structure_c.stableNorm());
    }
    diagnostics::neo_trace("candidate",
        "micro=%d,no=%ld,ns=%ld,hvp=%d,structure=%d,ro=%.17g,rs=%.17g,target=%.17g,shift=%.17g,radius=%.17g,step=%.17g,boundary=%d",
        result.iterations, static_cast<long>(qo.cols()),
        static_cast<long>(qs.cols()), result.orbital_actions,
        result.structure_actions, orbital_residual_norm, structure_residual_norm,
        kkt_residual_target, result.shift, options.trust_radius, result.step_norm,
        static_cast<int>(result.boundary));
    // Report already available Ritz data only. These are projected-model
    // eigenvalues, not a global spectral certificate or a new stopping test.
    if (diagnostics::neo_trace_enabled()) {
      diagnostics::neo_trace("spectrum",
          "micro=%d,hvp=%d,min=%.17g,max=%.17g,coordinate_residual=%.17g,physical_residual=%.17g",
          result.iterations, result.orbital_actions, spectrum_values_[0],
          spectrum_values_[spectrum_values_.size() - 1],
          result.kkt_residual.orbital.stableNorm(), orbital_residual_norm);
    }
    result.predicted_reduction =
        -problem_.orbital_gradient().dot(result.step.orbital) -
        0.5 * (result.step.orbital.dot(result.hessian_step.orbital) +
               result.step.structure.dot(result.hessian_step.structure));
    // Observation cannot accept a step or change solver tolerances. The
    // candidate is borrowed only for this call and is not yet certified.
    if (observe_candidate) observe_candidate(result);

    const bool stationary = result.residual_norm <= kkt_residual_target;
    const bool need_curvature_certificate =
        options.require_curvature_certificate || result.hard_case;
    Eigen::VectorXd curvature_orbital, curvature_structure;
    double curvature_orbital_target = 0.0;
    double curvature_structure_target = 0.0;
    bool curvature_converged = false;
    result.global_curvature_certified = false;
    // Optional Ritz diagnostics do not enter a regular KKT refinement. Defer
    // their five full-space image contractions, but evaluate them for every
    // final candidate and whenever the solver requires a curvature certificate.
    const auto measure_curvature = [&] {
      const Eigen::VectorXd minimum_coefficients = whitening *
          spectrum_vectors_.col(0);
      result.minimum_curvature_orbital.noalias() =
          qo * minimum_coefficients;
      const Eigen::VectorXd minimum_structure =
          response_coefficients * minimum_coefficients;
      const Eigen::VectorXd curvature_a = ao * minimum_coefficients;
      const Eigen::VectorXd curvature_bt = bt * minimum_structure;
      const Eigen::VectorXd curvature_metric =
          spectrum_values_[0] * mo * minimum_coefficients;
      const Eigen::VectorXd curvature_b = bo * minimum_coefficients;
      const Eigen::VectorXd curvature_c = cs * minimum_structure;
      curvature_orbital =
          curvature_a + curvature_bt - curvature_metric;
      curvature_structure = curvature_b + curvature_c;
      result.curvature_residual_norm = std::hypot(
          curvature_orbital.stableNorm(), curvature_structure.stableNorm());
      curvature_orbital_target =
          options.absolute_residual_tolerance +
          curvature_relative_tolerance * std::max({
              curvature_a.stableNorm(),
              curvature_bt.stableNorm(),
              curvature_metric.stableNorm()});
      curvature_structure_target =
          options.absolute_residual_tolerance +
          curvature_relative_tolerance * std::max(
              curvature_b.stableNorm(), curvature_c.stableNorm());
      result.curvature_residual_target = std::hypot(
          curvature_orbital_target, curvature_structure_target);
      curvature_converged =
          curvature_orbital.stableNorm() <= curvature_orbital_target &&
          curvature_structure.stableNorm() <= curvature_structure_target;
      result.global_curvature_certified =
          orbital_basis_.cols() == problem_.orbital_size() &&
          structure_basis_.cols() == problem_.structure_size() &&
          curvature_converged;
    };
    const int basis_size = static_cast<int>(
        orbital_basis_.cols() + structure_basis_.cols());
    if (need_curvature_certificate || stationary ||
        basis_size >= maximum_dimension) {
      measure_curvature();
    }
    const double positivity_tolerance =
        curvature_relative_tolerance * std::max({
            1.0, result.shift, std::abs(spectrum_values_[0])});
    const bool shifted_positive =
        spectrum_values_[0] + result.shift >= -positivity_tolerance;
    // A nonlinear keyframe certifies a regular boundary candidate by its
    // actual/predicted decrease.  Resolving an unrelated lowest Ritz vector
    // is therefore unnecessary unless explicitly requested.  A hard case is
    // different: its step itself depends on the minimum-curvature direction.
    const bool required_curvature_converged =
        !need_curvature_certificate || curvature_converged;
    const double descent = -problem_.orbital_gradient().dot(result.step.orbital);
    // Decide whether a nonlinear trial is ready before paying for an exact
    // residual response. TrialReady never claims physical KKT convergence;
    // the caller checks actual relaxed energy/gradient and disables this
    // branch near stationarity. Certified solves still perform the check below.
    const double coordinate_gradient = problem_.orbital_gradient().stableNorm();
    const double coordinate_target = options.absolute_residual_tolerance +
        coordinate_gradient * coordinate_gradient / (1.0 + coordinate_gradient);
    const double coordinate_residual = std::hypot(
        result.kkt_residual.orbital.stableNorm(), structure_residual_norm);
    if (allow_early_trial && !need_curvature_certificate && shifted_positive &&
        descent > 0.0 && basis_size > trial_basis_size_ &&
        coordinate_residual <= coordinate_target) {
      trial_basis_size_ = basis_size;
      result.stop_reason = NeoStopReason::TrialReady;
      diagnostics::neo_trace("historical_trial", "coordinate_residual=%.17g,coordinate_target=%.17g,physical_residual=%.17g,target=%.17g",
          coordinate_residual, coordinate_target, result.residual_norm, kkt_residual_target);
      return result;
    }
    // A certified solve must resolve the relaxed model, not first converge
    // its frozen-structure restriction. Seed the response of the actual
    // candidate once, reusing B p rather than solving for every B Q column.
    // Approximate nonlinear trials keep their deferred response policy.
    if (!allow_early_trial && !stationary && !response_seeded_ &&
        structure_basis_.cols() == 0 && basis_size < maximum_dimension &&
        problem_.has_structure_residual_solver() &&
        !result.kkt_residual.structure.isZero(0.0)) {
      response_seeded_ = true;
      const Eigen::VectorXd response = problem_.solve_structure_residual(
          result.kkt_residual.structure);
      if (!response.isZero(0.0)) {
        const auto image = problem_.apply_structure_coupling(response);
        ++result.coupled_actions;
        ++result.structure_actions;
        if (retain_structure_image(response, image)) continue;
      }
    }
    if (stationary && required_curvature_converged && shifted_positive) {
      // A small raw structure residual can still produce a large orbital
      // response when C is ill-conditioned. Check the actual response before
      // certifying the projected candidate.
      if (problem_.has_structure_residual_solver() &&
          !result.kkt_residual.structure.isZero(0.0)) {
        const Eigen::VectorXd correction = problem_.solve_structure_residual(
            result.kkt_residual.structure);
        if (!correction.isZero(0.0)) {
          const ResponseNeoDirection image =
              problem_.apply_structure_coupling(correction);
          ++result.coupled_actions;
          ++result.structure_actions;
          const Eigen::VectorXd corrected_orbital =
              result.kkt_residual.orbital + image.orbital;
          const Eigen::VectorXd corrected_structure =
              result.kkt_residual.structure + image.structure;
          const double corrected_norm = std::hypot(
              problem_.orbital_dual_norm(corrected_orbital),
              corrected_structure.stableNorm());
          diagnostics::neo_trace("response_check", "before=%.17g,after=%.17g,target=%.17g,retained=%d",
              result.residual_norm, corrected_norm, kkt_residual_target,
              static_cast<int>(corrected_norm > kkt_residual_target));
          if (corrected_norm > kkt_residual_target) {
            const int basis_size = static_cast<int>(
                orbital_basis_.cols() + structure_basis_.cols());
            if (structure_basis_.cols() == problem_.structure_size()) {
              result.stop_reason = NeoStopReason::NumericalFailure;
              return result;
            }
            if (basis_size >= maximum_dimension) {
              result.stop_reason = NeoStopReason::SubspaceLimit;
              return result;
            }
            if (!retain_structure_image(correction, image)) {
              result.stop_reason = NeoStopReason::NumericalFailure;
              return result;
            }
            continue;
          }
          result.step.structure += correction;
          result.hessian_step.orbital += image.orbital;
          result.hessian_step.structure += image.structure;
          result.kkt_residual.orbital = corrected_orbital;
          result.kkt_residual.structure = corrected_structure;
          result.residual_norm = corrected_norm;
          result.predicted_reduction =
              -problem_.orbital_gradient().dot(result.step.orbital) -
              0.5 * (result.step.orbital.dot(result.hessian_step.orbital) +
                     result.step.structure.dot(result.hessian_step.structure));
        }
      }
      result.stop_reason = NeoStopReason::Converged;
      return result;
    }

    // CIAH's early-start test applies to the normalized augmented vector,
    // not the unscaled Newton residual. For [1, alpha*p, alpha*q], the top
    // AH equation gives alpha^2 = lambda / (-g^T p). The orbital trust
    // metric yields normalization sqrt(1 + alpha^2 p^T G p); the structure
    // coordinates have no trust penalty. This candidate is NOT certified.
    // Three directions and the 2.5 start threshold follow PySCF mc1step;
    // exact energy/gradient keyframes decide whether it is usable.
    if (allow_early_trial && !need_curvature_certificate && shifted_positive &&
        result.boundary && result.shift > 0.0 && descent > 0.0 &&
        qo.cols() >= 3 && basis_size > trial_basis_size_) {
      const double augmented_scale = std::sqrt(result.shift / descent);
      const double augmented_residual = result.residual_norm * augmented_scale /
          std::hypot(1.0, augmented_scale * result.step_norm);
      if (augmented_residual < 2.5) {
        trial_basis_size_ = basis_size;
        result.stop_reason = NeoStopReason::TrialReady;
        diagnostics::neo_trace("early_trial", "residual=%.17g,augmented_residual=%.17g,target=%.17g",
            result.residual_norm, augmented_residual, kkt_residual_target);
        return result;
      }
    }

    if (basis_size >= maximum_dimension) {
      result.stop_reason = NeoStopReason::SubspaceLimit;
      return result;
    }

    // Refine only the structure response required by the current candidate.
    residual_phase.stop();
    projected_phase.stop();
    // Closing C Z = -B Q for every retained orbital basis column performs
    // response solves that do not enter p = Q y.  The complete coupled KKT
    // residual below certifies the direction-local response without that
    // basis-wide work.
    bool expanded = false;
    diagnostics::neo_trace_context.source = "structure_residual";
    if (!stationary && structure_residual_norm >= orbital_residual_norm) {
      const bool add_curvature = need_curvature_certificate &&
          curvature_structure.stableNorm() > curvature_structure_target;
      Eigen::MatrixXd directions(
          problem_.structure_size(), add_curvature ? 2 : 1);
      directions.col(0) = -result.kkt_residual.structure;
      if (add_curvature) directions.col(1) = -curvature_structure;
      expanded = append_structure_block(
          std::move(directions), &result.coupled_actions,
          &result.structure_actions);
    }
    if (!expanded && need_curvature_certificate &&
        curvature_structure.stableNorm() >
        curvature_structure_target) {
      diagnostics::neo_trace_context.source = "structure_curvature";
      expanded = append_structure(
          -curvature_structure, &result.coupled_actions,
          &result.structure_actions);
    }
    if (!expanded && !stationary) {
      diagnostics::neo_trace_context.source = "orbital_residual";
      expanded = append_orbital(
          preconditioned_orbital(
              problem_, -result.kkt_residual.orbital,
              relative_preconditioner_budget(
                  kkt_residual_target, orbital_residual_norm),
              result.shift),
          &result.coupled_actions,
          &result.orbital_actions);
    }
    if (!expanded && !stationary && structure_residual_norm > 0.0) {
      diagnostics::neo_trace_context.source = "structure_residual";
      expanded = append_structure(
          -result.kkt_residual.structure, &result.coupled_actions,
          &result.structure_actions);
    }
    if (!expanded && need_curvature_certificate &&
        curvature_orbital.stableNorm() >
        curvature_orbital_target) {
      diagnostics::neo_trace_context.source = "orbital_curvature";
      expanded = append_orbital(
          preconditioned_orbital(
              problem_, -curvature_orbital,
              relative_preconditioner_budget(
                  kkt_residual_target,
                  problem_.orbital_dual_norm(curvature_orbital))),
          &result.coupled_actions,
          &result.orbital_actions);
    }
    while (!expanded && orbital_canonical_ < problem_.orbital_size()) {
      diagnostics::neo_trace_context.source = "orbital_canonical";
      Eigen::VectorXd direction =
          Eigen::VectorXd::Zero(problem_.orbital_size());
      direction[orbital_canonical_++] = 1.0;
      expanded = append_orbital(
          std::move(direction), &result.coupled_actions,
          &result.orbital_actions);
    }
    while (!expanded && structure_canonical_ < problem_.structure_size()) {
      diagnostics::neo_trace_context.source = "structure_canonical";
      Eigen::VectorXd direction =
          Eigen::VectorXd::Zero(problem_.structure_size());
      direction[structure_canonical_++] = 1.0;
      expanded = append_structure(
          std::move(direction), &result.coupled_actions,
          &result.structure_actions);
    }
    diagnostics::neo_trace("expand", "source=%s,success=%d,no=%ld,ns=%ld",
        diagnostics::neo_trace_context.source, static_cast<int>(expanded),
        static_cast<long>(orbital_basis_.cols()),
        static_cast<long>(structure_basis_.cols()));
    if (!expanded) {
      if (!need_curvature_certificate && !stationary) measure_curvature();
      result.stop_reason = NeoStopReason::NumericalFailure;
      return result;
    }
  }
}

ResponseNeoResult solve_response_neo(
    const ResponseNeoProblem& problem,
    const NeoOptions& options) {
  ResponseNeoWorkspace workspace(problem);
  return workspace.solve(options);
}

}  // namespace xmvb::vb
