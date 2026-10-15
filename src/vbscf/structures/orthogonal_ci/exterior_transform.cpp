#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iterator>
#include <list>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <unordered_map>

#include "core/openmp.hpp"
#include "core/timing.hpp"
#include "mcscf/ci/product_support.hpp"

namespace xmvb::vb {
namespace {

std::size_t binomial(int n, int k) {
  if (k < 0 || k > n) {
    return 0;
  }
  k = std::min(k, n - k);
  std::size_t value = 1;
  for (int index = 1; index <= k; ++index) {
    value = value * static_cast<std::size_t>(n - k + index) /
        static_cast<std::size_t>(index);
  }
  return value;
}

double replacement_sign(std::uint64_t source, int removed, int inserted) {
  const int low = std::min(removed, inserted);
  const int high = std::max(removed, inserted);
  std::uint64_t between = 0;
  if (high - low > 1) {
    between = source &
        (((std::uint64_t{1} << high) - 1) ^
         ((std::uint64_t{1} << (low + 1)) - 1));
  }
  return __builtin_popcountll(between) % 2 == 0 ? 1.0 : -1.0;
}

}  // namespace

struct ExteriorOrbitalTransform::RowCache {
  struct Entry {
    std::uint64_t key = 0;
    std::shared_ptr<const Eigen::MatrixXd> rows;
    std::size_t bytes = 0;
  };

  std::mutex mutex;
  std::list<Entry> entries;
  std::unordered_map<std::uint64_t, std::list<Entry>::iterator> index;
  std::atomic<std::size_t> cached_bytes{0};
};

ExteriorOrbitalTransform::ExteriorOrbitalTransform(
    const std::vector<std::vector<int>>& determinants,
    const Eigen::Ref<const Eigen::MatrixXd>& upper_orbital_transform)
    : n_orbitals_(static_cast<int>(upper_orbital_transform.rows())),
      orbital_transform_(upper_orbital_transform),
      row_cache_(std::make_shared<RowCache>()) {
  if (n_orbitals_ <= 0 || n_orbitals_ > 63 ||
      upper_orbital_transform.cols() != n_orbitals_ ||
      determinants.empty() || !upper_orbital_transform.allFinite()) {
    throw std::invalid_argument("invalid exterior orbital transform inputs");
  }
  const double lower_error =
      upper_orbital_transform.template triangularView<Eigen::StrictlyLower>()
          .toDenseMatrix()
          .cwiseAbs()
          .maxCoeff();
  const double scale = std::max(
      1.0,
      upper_orbital_transform.cwiseAbs().maxCoeff());
  if (lower_error > 64.0 * std::numeric_limits<double>::epsilon() * scale ||
      (upper_orbital_transform.diagonal().array() <= 0.0).any()) {
    throw std::invalid_argument(
        "exterior orbital transform must be nonsingular upper triangular");
  }

  n_electrons_ = static_cast<int>(determinants.front().size());
  if (determinants.size() != binomial(n_orbitals_, n_electrons_)) {
    throw std::invalid_argument(
        "exterior transform requires a complete fixed-spin space");
  }
  masks_.reserve(determinants.size());
  std::unordered_map<std::uint64_t, int> index_by_mask;
  index_by_mask.reserve(determinants.size());
  for (int index = 0; index < static_cast<int>(determinants.size()); ++index) {
    const auto& occupied = determinants[index];
    if (static_cast<int>(occupied.size()) != n_electrons_) {
      throw std::invalid_argument(
          "fixed-spin determinants have inconsistent electron counts");
    }
    std::uint64_t mask = 0;
    int previous = -1;
    for (const int orbital : occupied) {
      if (orbital <= previous || orbital < 0 || orbital >= n_orbitals_) {
        throw std::invalid_argument(
            "fixed-spin determinant occupation is invalid");
      }
      mask |= std::uint64_t{1} << orbital;
      previous = orbital;
    }
    if (!index_by_mask.emplace(mask, index).second) {
      throw std::invalid_argument("fixed-spin determinant is duplicated");
    }
    masks_.push_back(mask);
  }

  determinant_scales_.resize(static_cast<int>(masks_.size()));
  for (int determinant = 0;
       determinant < static_cast<int>(masks_.size());
       ++determinant) {
    double scale_factor = 1.0;
    for (int orbital = 0; orbital < n_orbitals_; ++orbital) {
      if ((masks_[determinant] & (std::uint64_t{1} << orbital)) != 0) {
        scale_factor *= upper_orbital_transform(orbital, orbital);
      }
    }
    determinant_scales_[determinant] = scale_factor;
  }

  // R = U D, where D contains the column diagonal and U is unit upper
  // triangular. Reduce U to identity by right-looking column eliminations.
  // Applying the inverses in elimination order evaluates U on a vector
  // because matrix factors act from right to left.
  Eigen::MatrixXd unit_upper = upper_orbital_transform;
  for (int column = 0; column < n_orbitals_; ++column) {
    unit_upper.col(column) /= upper_orbital_transform(column, column);
  }
  struct OrbitalShear {
    int inserted = 0;
    int removed = 0;
    double coefficient = 0.0;
  };
  std::vector<OrbitalShear> eliminations;
  for (int removed = 1; removed < n_orbitals_; ++removed) {
    for (int inserted = removed - 1; inserted >= 0; --inserted) {
      const double coefficient = unit_upper(inserted, removed);
      eliminations.push_back({inserted, removed, coefficient});
      unit_upper.col(removed).noalias() -=
          coefficient * unit_upper.col(inserted);
    }
  }

  shears_.reserve(eliminations.size());
  for (const OrbitalShear& operation : eliminations) {
    Shear shear;
    shear.coefficient = operation.coefficient;
    const std::uint64_t removed_bit =
        std::uint64_t{1} << operation.removed;
    const std::uint64_t inserted_bit =
        std::uint64_t{1} << operation.inserted;
    for (int source = 0; source < static_cast<int>(masks_.size()); ++source) {
      const std::uint64_t source_mask = masks_[source];
      if ((source_mask & removed_bit) == 0 ||
          (source_mask & inserted_bit) != 0) {
        continue;
      }
      const std::uint64_t target_mask =
          (source_mask ^ removed_bit) | inserted_bit;
      const auto target = index_by_mask.find(target_mask);
      if (target == index_by_mask.end()) {
        throw std::invalid_argument(
            "exterior transform requires a complete fixed-spin space");
      }
      shear.pairs.push_back(DeterminantPair{
          source,
          target->second,
          replacement_sign(
              source_mask,
              operation.removed,
              operation.inserted)});
    }
    shears_.push_back(std::move(shear));
  }
}

void ExteriorOrbitalTransform::validate_left(
    const Eigen::MatrixXd& coefficients) const {
  if (coefficients.rows() != dimension()) {
    throw std::invalid_argument(
        "left exterior transform determinant dimension mismatch");
  }
}

void ExteriorOrbitalTransform::validate_right(
    const Eigen::MatrixXd& coefficients) const {
  if (coefficients.cols() != dimension()) {
    throw std::invalid_argument(
        "right exterior transform determinant dimension mismatch");
  }
}

void ExteriorOrbitalTransform::apply_left(
    Eigen::MatrixXd* coefficients) const {
  if (coefficients == nullptr) {
    throw std::invalid_argument("left exterior transform output is null");
  }
  validate_left(*coefficients);
  transform_left(*coefficients, false, false);
}

void ExteriorOrbitalTransform::apply_left_block(Eigen::Ref<Eigen::MatrixXd> coefficients) const {
  if (coefficients.rows() != dimension())
    throw std::invalid_argument("left exterior block dimension mismatch");
  transform_left(coefficients, false, false);
}

void ExteriorOrbitalTransform::apply_right(
    Eigen::MatrixXd* coefficients) const {
  if (coefficients == nullptr) {
    throw std::invalid_argument("right exterior transform output is null");
  }
  apply_right_block(*coefficients);
}

void ExteriorOrbitalTransform::apply_right_block(
    Eigen::Ref<Eigen::MatrixXd> coefficients) const {
  if (coefficients.cols() != dimension()) {
    throw std::invalid_argument(
        "right exterior transform determinant dimension mismatch");
  }
  transform_right(coefficients, false, false);
}

void ExteriorOrbitalTransform::apply_adjoint_left(
    Eigen::MatrixXd* coefficients) const {
  if (coefficients == nullptr) {
    throw std::invalid_argument("left exterior adjoint output is null");
  }
  validate_left(*coefficients);
  transform_left(*coefficients, true, false);
}

ExteriorProjection ExteriorOrbitalTransform::projection(const std::vector<int>& rows) const {
  std::vector<unsigned char> needed(dimension(), 0);
  for (int row : rows) {
    if (row < 0 || row >= dimension()) throw std::invalid_argument("exterior projection row is outside the spin space");
    needed[row] = 1;
  }
  ExteriorProjection result;
  result.shear_pairs.resize(shears_.size());
  // The adjoint executes shears in reverse. Walk their dependency graph in
  // forward order, starting from the requested final rows. Pairs in one
  // elementary shear are disjoint, so this sweep needs no transitive iteration.
  for (std::size_t i = 0; i < shears_.size(); ++i) {
    const auto& pairs = shears_[i].pairs;
    for (int j = 0; j < static_cast<int>(pairs.size()); ++j)
      if (needed[pairs[j].source]) {
        result.shear_pairs[i].push_back(j);
        needed[pairs[j].target] = 1;
      }
  }
  return result;
}

void ExteriorOrbitalTransform::apply_adjoint_left(
    Eigen::MatrixXd* coefficients, const ExteriorProjection& projection) const {
  if (!coefficients)
    throw std::invalid_argument("invalid exterior projection plan");
  apply_adjoint_left_block(*coefficients, projection);
}

void ExteriorOrbitalTransform::apply_adjoint_left_block(
    Eigen::Ref<Eigen::MatrixXd> coefficients, const ExteriorProjection& projection) const {
  if (coefficients.rows() != dimension() || projection.shear_pairs.size() != shears_.size())
    throw std::invalid_argument("invalid exterior projection block");
  profiling::ScopedPhase phase(profiling::Phase::CiTransformLeft);
  const int workers = bounded_openmp_thread_count(
      coarse_openmp_thread_count(effective_openmp_thread_count(), coefficients.size()), coefficients.cols());
#pragma omp parallel for schedule(static) if(workers > 1) num_threads(workers)
  for (Eigen::Index column = 0; column < coefficients.cols(); ++column) {
    auto values = coefficients.col(column);
    for (int i = static_cast<int>(shears_.size()) - 1; i >= 0; --i) {
      const auto& shear = shears_[i];
      for (int j : projection.shear_pairs[i]) {
        const auto& pair = shear.pairs[j];
        values[pair.source] += shear.coefficient * pair.sign * values[pair.target];
      }
    }
    values.array() *= determinant_scales_.array();
  }
}

void ExteriorOrbitalTransform::apply_supported_left(
    Eigen::MatrixXd* coefficients, const mcscf::ProductSupport& support,
    int beta_begin, int beta_count, bool adjoint) const {
  xmvb::profiling::ScopedPhase phase(xmvb::profiling::Phase::CiTransformLeft);
  validate_left(*coefficients);
  const int workers = bounded_openmp_thread_count(
      coarse_openmp_thread_count(effective_openmp_thread_count(), coefficients->size()),
      coefficients->cols());
  const int n_shears = static_cast<int>(shears_.size());
#pragma omp parallel if(workers > 1) num_threads(workers)
  {
    std::vector<std::uint8_t> allowed(dimension());
#pragma omp for schedule(static)
    for (Eigen::Index column = 0; column < coefficients->cols(); ++column) {
      auto values = coefficients->col(column);
      const int beta = beta_begin + column % beta_count;
      std::fill(allowed.begin(), allowed.end(), 0);
      auto paths = support.column(beta, 0);
      int alpha_index, beta_index;
      bool in_support;
      while (paths.next(alpha_index, beta_index, in_support)) allowed[alpha_index] = 1;
      for (int alpha = 0; alpha < dimension(); ++alpha) {
        if (!allowed[alpha]) values[alpha] = 0.0;
        else if (!adjoint) values[alpha] *= determinant_scales_[alpha];
      }
      for (int index = 0; index < n_shears; ++index) {
        const Shear& shear = shears_[adjoint ? n_shears - 1 - index : index];
        for (const auto& pair : shear.pairs) {
          if (!allowed[pair.source] || !allowed[pair.target]) continue;
          const int target = adjoint ? pair.source : pair.target;
          const int source = adjoint ? pair.target : pair.source;
          values[target] += shear.coefficient * pair.sign * values[source];
        }
      }
      if (adjoint)
        for (int alpha = 0; alpha < dimension(); ++alpha)
          if (allowed[alpha]) values[alpha] *= determinant_scales_[alpha];
    }
  }
}

void ExteriorOrbitalTransform::apply_adjoint_right(
    Eigen::MatrixXd* coefficients) const {
  if (coefficients == nullptr) {
    throw std::invalid_argument("right exterior adjoint output is null");
  }
  apply_adjoint_right_block(*coefficients);
}

void ExteriorOrbitalTransform::apply_adjoint_right_block(
    Eigen::Ref<Eigen::MatrixXd> coefficients) const {
  if (coefficients.cols() != dimension()) {
    throw std::invalid_argument(
        "right exterior transform determinant dimension mismatch");
  }
  transform_right(coefficients, true, false);
}

std::shared_ptr<const Eigen::MatrixXd> ExteriorOrbitalTransform::row_block(
    int begin, int count) const {
  if (begin < 0 || count <= 0 || begin > dimension() - count) {
    throw std::invalid_argument("exterior row block is out of range");
  }
  constexpr std::size_t kCacheBytes = 64ULL * 1024ULL * 1024ULL;
  const std::uint64_t key =
      (static_cast<std::uint64_t>(static_cast<std::uint32_t>(begin)) << 32) |
      static_cast<std::uint32_t>(count);
  {
    std::lock_guard<std::mutex> lock(row_cache_->mutex);
    const auto found = row_cache_->index.find(key);
    if (found != row_cache_->index.end()) {
      row_cache_->entries.splice(row_cache_->entries.begin(),
                                 row_cache_->entries, found->second);
      return found->second->rows;
    }
  }

  // Generate outside the mutex; other readers may use cached rows meanwhile.
  auto rows = std::make_shared<Eigen::MatrixXd>(count, dimension());
  rows->setZero();
  rows->middleCols(begin, count).setIdentity();
  apply_adjoint_right_block(*rows);
  const std::size_t width = static_cast<std::size_t>(dimension());
  if (static_cast<std::size_t>(count) >
      std::numeric_limits<std::size_t>::max() / width / sizeof(double)) {
    throw std::length_error("exterior row block byte count overflows");
  }
  const std::size_t bytes =
      static_cast<std::size_t>(count) * width * sizeof(double);
  if (bytes > kCacheBytes) return rows;

  std::lock_guard<std::mutex> lock(row_cache_->mutex);
  const auto found = row_cache_->index.find(key);
  if (found != row_cache_->index.end()) {
    row_cache_->entries.splice(row_cache_->entries.begin(),
                               row_cache_->entries, found->second);
    return found->second->rows;
  }
  std::size_t resident = row_cache_->cached_bytes.load(std::memory_order_relaxed);
  while (resident > kCacheBytes - bytes) {
    auto last = std::prev(row_cache_->entries.end());
    resident -= last->bytes;
    row_cache_->index.erase(last->key);
    row_cache_->entries.erase(last);
  }
  row_cache_->cached_bytes.store(resident, std::memory_order_relaxed);
  row_cache_->entries.push_front(RowCache::Entry{key, rows, bytes});
  try {
    row_cache_->index.emplace(key, row_cache_->entries.begin());
  } catch (...) {
    row_cache_->entries.pop_front();
    throw;
  }
  row_cache_->cached_bytes.store(resident + bytes, std::memory_order_relaxed);
  return rows;
}

void ExteriorOrbitalTransform::validate_inverse() const {
  if (!determinant_scales_.allFinite() ||
      (determinant_scales_.array() == 0.0).any()) {
    throw std::runtime_error("exterior transform has noninvertible numerical scales");
  }
}

void ExteriorOrbitalTransform::apply_inverse_left(
    Eigen::MatrixXd* coefficients) const {
  if (coefficients == nullptr)
    throw std::invalid_argument("inverse exterior transform output is null");
  validate_left(*coefficients);
  validate_inverse();
  transform_left(*coefficients, false, true);
}

void ExteriorOrbitalTransform::apply_inverse_adjoint_left(
    Eigen::MatrixXd* coefficients) const {
  if (coefficients == nullptr)
    throw std::invalid_argument("inverse exterior adjoint output is null");
  validate_left(*coefficients);
  validate_inverse();
  transform_left(*coefficients, true, true);
}

void ExteriorOrbitalTransform::apply_inverse_right_block(
    Eigen::Ref<Eigen::MatrixXd> coefficients) const {
  if (coefficients.cols() != dimension())
    throw std::invalid_argument("inverse exterior column dimension mismatch");
  validate_inverse();
  transform_right(coefficients, false, true);
}

void ExteriorOrbitalTransform::apply_inverse_adjoint_right_block(
    Eigen::Ref<Eigen::MatrixXd> coefficients) const {
  if (coefficients.cols() != dimension())
    throw std::invalid_argument("inverse exterior adjoint column dimension mismatch");
  validate_inverse();
  transform_right(coefficients, true, true);
}

void ExteriorOrbitalTransform::transform_left(
    Eigen::Ref<Eigen::MatrixXd> coefficients, bool adjoint, bool inverse) const {
  xmvb::profiling::ScopedPhase phase(xmvb::profiling::Phase::CiTransformLeft);
  const bool reverse = adjoint != inverse;
  const int n_shears = static_cast<int>(shears_.size());
  const int workers = bounded_openmp_thread_count(
      coarse_openmp_thread_count(effective_openmp_thread_count(), coefficients.size()),
      coefficients.cols());
  // A column is an independent vector. Keep it resident while applying the
  // complete ordered product, rather than streaming strided rows per shear.
#pragma omp parallel for schedule(static) if(workers > 1) num_threads(workers)
  for (Eigen::Index column = 0; column < coefficients.cols(); ++column) {
    auto values = coefficients.col(column);
    if (!reverse) {
      if (inverse) values.array() /= determinant_scales_.array();
      else values.array() *= determinant_scales_.array();
    }
    for (int index = 0; index < n_shears; ++index) {
      const Shear& shear = shears_[reverse ? n_shears - 1 - index : index];
      const double factor = inverse ? -shear.coefficient : shear.coefficient;
      for (const DeterminantPair& pair : shear.pairs) {
        const int target = adjoint ? pair.source : pair.target;
        const int source = adjoint ? pair.target : pair.source;
        values[target] += (factor * pair.sign) * values[source];
      }
    }
    if (reverse) {
      if (inverse) values.array() /= determinant_scales_.array();
      else values.array() *= determinant_scales_.array();
    }
  }
}

void ExteriorOrbitalTransform::transform_right(
    Eigen::Ref<Eigen::MatrixXd> coefficients, bool adjoint, bool inverse) const {
  xmvb::profiling::ScopedPhase phase(xmvb::profiling::Phase::CiTransformRight);
  const bool reverse = adjoint != inverse;
  const int n_shears = static_cast<int>(shears_.size());
  const int workers = bounded_openmp_thread_count(
      coarse_openmp_thread_count(effective_openmp_thread_count(), coefficients.size()),
      coefficients.rows());
  // Workers own disjoint row ranges. Each range executes all shears without
  // intermediate barriers. Column segments stay contiguous for SIMD updates.
#pragma omp parallel for schedule(static) if(workers > 1) num_threads(workers)
  for (int part = 0; part < workers; ++part) {
    const Eigen::Index first = coefficients.rows() * part / workers;
    const Eigen::Index end = coefficients.rows() * (part + 1) / workers;
    auto values = coefficients.middleRows(first, end - first);
    if (!reverse) {
      if (inverse) values.array().rowwise() /= determinant_scales_.transpose().array();
      else values.array().rowwise() *= determinant_scales_.transpose().array();
    }
    for (int index = 0; index < n_shears; ++index) {
      const Shear& shear = shears_[reverse ? n_shears - 1 - index : index];
      const double factor = inverse ? -shear.coefficient : shear.coefficient;
      for (const DeterminantPair& pair : shear.pairs) {
        const int target = adjoint ? pair.source : pair.target;
        const int source = adjoint ? pair.target : pair.source;
        values.col(target) += (factor * pair.sign) * values.col(source);
      }
    }
    if (reverse) {
      if (inverse) values.array().rowwise() /= determinant_scales_.transpose().array();
      else values.array().rowwise() *= determinant_scales_.transpose().array();
    }
  }
}

ExteriorTransformDirection ExteriorOrbitalTransform::direction(
    const Eigen::Ref<const Eigen::MatrixXd>& delta_orbital_transform) const {
  if (delta_orbital_transform.rows() != n_orbitals_ ||
      delta_orbital_transform.cols() != n_orbitals_ ||
      !delta_orbital_transform.allFinite()) {
    throw std::invalid_argument(
        "exterior transform direction has incompatible dimensions");
  }
  const double lower_error =
      delta_orbital_transform
          .template triangularView<Eigen::StrictlyLower>()
          .toDenseMatrix()
          .cwiseAbs()
          .maxCoeff();
  const double scale = std::max(
      1.0,
      delta_orbital_transform.cwiseAbs().maxCoeff());
  if (lower_error > 64.0 * std::numeric_limits<double>::epsilon() * scale) {
    throw std::invalid_argument(
        "exterior transform direction must be upper triangular");
  }

  ExteriorTransformDirection result;
  result.determinant_scales.resize(dimension());
  for (int determinant = 0; determinant < dimension(); ++determinant) {
    double logarithmic_direction = 0.0;
    for (int orbital = 0; orbital < n_orbitals_; ++orbital) {
      if ((masks_[determinant] & (std::uint64_t{1} << orbital)) != 0) {
        logarithmic_direction +=
            delta_orbital_transform(orbital, orbital) /
            orbital_transform_(orbital, orbital);
      }
    }
    result.determinant_scales[determinant] =
        determinant_scales_[determinant] * logarithmic_direction;
  }

  Eigen::MatrixXd unit_upper = orbital_transform_;
  Eigen::MatrixXd delta_unit_upper = delta_orbital_transform;
  for (int column = 0; column < n_orbitals_; ++column) {
    const double diagonal = orbital_transform_(column, column);
    const double delta_diagonal = delta_orbital_transform(column, column);
    delta_unit_upper.col(column) =
        delta_unit_upper.col(column) / diagonal -
        unit_upper.col(column) * (delta_diagonal / (diagonal * diagonal));
    unit_upper.col(column) /= diagonal;
  }

  result.shear_coefficients.reserve(shears_.size());
  for (int removed = 1; removed < n_orbitals_; ++removed) {
    for (int inserted = removed - 1; inserted >= 0; --inserted) {
      const double coefficient = unit_upper(inserted, removed);
      const double delta_coefficient =
          delta_unit_upper(inserted, removed);
      result.shear_coefficients.push_back(delta_coefficient);
      delta_unit_upper.col(removed).noalias() -=
          delta_coefficient * unit_upper.col(inserted) +
          coefficient * delta_unit_upper.col(inserted);
      unit_upper.col(removed).noalias() -=
          coefficient * unit_upper.col(inserted);
    }
  }
  if (result.shear_coefficients.size() != shears_.size()) {
    throw std::logic_error(
        "exterior transform direction does not match its shear schedule");
  }
  return result;
}

void ExteriorOrbitalTransform::apply_directional_left(
    const ExteriorTransformDirection& direction,
    Eigen::MatrixXd* coefficients,
    Eigen::MatrixXd* directional_coefficients) const {
  transform_directional_left(direction, coefficients, directional_coefficients, false);
}

void ExteriorOrbitalTransform::transform_directional_left(
    const ExteriorTransformDirection& direction,
    Eigen::MatrixXd* coefficients, Eigen::MatrixXd* tangent,
    bool adjoint) const {
  xmvb::profiling::ScopedPhase phase(xmvb::profiling::Phase::CiDirectionalTransform);
  if (coefficients == nullptr || tangent == nullptr || coefficients == tangent ||
      coefficients->rows() != tangent->rows() ||
      coefficients->cols() != tangent->cols() ||
      direction.determinant_scales.size() != dimension() ||
      direction.shear_coefficients.size() != shears_.size()) {
    throw std::invalid_argument("invalid directional left exterior action");
  }
  validate_left(*coefficients);
  const int n_shears = static_cast<int>(shears_.size());
  const int workers = bounded_openmp_thread_count(
      coarse_openmp_thread_count(effective_openmp_thread_count(), coefficients->size()),
      coefficients->cols());
  // Different columns are independent. Primal and tangent follow the same
  // ordered shears in one pass, with no thread-local copies or shear barriers.
#pragma omp parallel for schedule(static) if(workers > 1) num_threads(workers)
  for (Eigen::Index column = 0; column < coefficients->cols(); ++column) {
    auto value = coefficients->col(column);
    auto delta = tangent->col(column);
    const auto scale = [&]() {
      delta.array() = determinant_scales_.array() * delta.array() +
                     direction.determinant_scales.array() * value.array();
      value.array() *= determinant_scales_.array();
    };
    if (!adjoint) scale();
    for (int index = 0; index < n_shears; ++index) {
      const int ordered = adjoint ? n_shears - 1 - index : index;
      const Shear& shear = shears_[ordered];
      const double dc = direction.shear_coefficients[ordered];
      for (const DeterminantPair& pair : shear.pairs) {
        const int source = adjoint ? pair.target : pair.source;
        const int target = adjoint ? pair.source : pair.target;
        delta[target] += pair.sign *
            (dc * value[source] + shear.coefficient * delta[source]);
        value[target] += (shear.coefficient * pair.sign) * value[source];
      }
    }
    if (adjoint) scale();
  }
}

void ExteriorOrbitalTransform::apply_directional_right(
    const ExteriorTransformDirection& direction,
    Eigen::MatrixXd* coefficients,
    Eigen::MatrixXd* directional_coefficients) const {
  xmvb::profiling::ScopedPhase phase(xmvb::profiling::Phase::CiDirectionalTransform);
  if (coefficients == nullptr || directional_coefficients == nullptr ||
      coefficients->rows() != directional_coefficients->rows() ||
      coefficients->cols() != directional_coefficients->cols() ||
      direction.determinant_scales.size() != dimension() ||
      direction.shear_coefficients.size() != shears_.size()) {
    throw std::invalid_argument("invalid directional right exterior action");
  }
  validate_right(*coefficients);
  for (int column = 0; column < dimension(); ++column) {
    directional_coefficients->col(column) =
        determinant_scales_[column] * directional_coefficients->col(column) +
        direction.determinant_scales[column] * coefficients->col(column);
    coefficients->col(column) *= determinant_scales_[column];
  }
  for (std::size_t shear_index = 0;
       shear_index < shears_.size();
       ++shear_index) {
    const Shear& shear = shears_[shear_index];
    const double delta_coefficient =
        direction.shear_coefficients[shear_index];
    for (const DeterminantPair& pair : shear.pairs) {
      directional_coefficients->col(pair.target) += pair.sign *
          (delta_coefficient * coefficients->col(pair.source) +
           shear.coefficient * directional_coefficients->col(pair.source));
      coefficients->col(pair.target) +=
          (shear.coefficient * pair.sign) * coefficients->col(pair.source);
    }
  }
}

void ExteriorOrbitalTransform::apply_directional_adjoint_left(
    const ExteriorTransformDirection& direction,
    Eigen::MatrixXd* coefficients,
    Eigen::MatrixXd* directional_coefficients) const {
  transform_directional_left(direction, coefficients, directional_coefficients, true);
}

void ExteriorOrbitalTransform::apply_directional_adjoint_right(
    const ExteriorTransformDirection& direction,
    Eigen::MatrixXd* coefficients,
    Eigen::MatrixXd* directional_coefficients) const {
  xmvb::profiling::ScopedPhase phase(xmvb::profiling::Phase::CiDirectionalTransform);
  if (coefficients == nullptr || directional_coefficients == nullptr ||
      coefficients->rows() != directional_coefficients->rows() ||
      coefficients->cols() != directional_coefficients->cols() ||
      direction.determinant_scales.size() != dimension() ||
      direction.shear_coefficients.size() != shears_.size()) {
    throw std::invalid_argument("invalid directional adjoint-right exterior action");
  }
  validate_right(*coefficients);
  for (std::size_t reverse = shears_.size(); reverse-- > 0;) {
    const Shear& shear = shears_[reverse];
    const double delta_coefficient = direction.shear_coefficients[reverse];
    for (const DeterminantPair& pair : shear.pairs) {
      directional_coefficients->col(pair.source) += pair.sign *
          (delta_coefficient * coefficients->col(pair.target) +
           shear.coefficient * directional_coefficients->col(pair.target));
      coefficients->col(pair.source) +=
          (shear.coefficient * pair.sign) * coefficients->col(pair.target);
    }
  }
  for (int column = 0; column < dimension(); ++column) {
    directional_coefficients->col(column) =
        determinant_scales_[column] * directional_coefficients->col(column) +
        direction.determinant_scales[column] * coefficients->col(column);
    coefficients->col(column) *= determinant_scales_[column];
  }
}

std::size_t ExteriorOrbitalTransform::shear_pair_count() const noexcept {
  std::size_t count = 0;
  for (const Shear& shear : shears_) {
    count += shear.pairs.size();
  }
  return count;
}

Eigen::MatrixXd ExteriorOrbitalTransform::generator_pullback(
    const Eigen::MatrixXd& seed, const Eigen::MatrixXd& image,
    const Eigen::MatrixXd* delta_seed, const Eigen::MatrixXd* delta_image,
    const Eigen::MatrixXd* delta_upper) const {
  const bool directional = delta_upper != nullptr;
  if (seed.rows() != image.rows() || seed.cols() != image.cols() ||
      image.rows() != dimension() ||
      directional != (delta_seed != nullptr) || directional != (delta_image != nullptr))
    throw std::invalid_argument("invalid generator pullback panels");
  Eigen::MatrixXd g = Eigen::MatrixXd::Zero(n_orbitals_, n_orbitals_);
  Eigen::MatrixXd dg = g;
  const auto contract = [&](int target, int source, bool tangent) -> double {
    if (tangent) return delta_seed->row(target).dot(image.row(source)) +
        seed.row(target).dot(delta_image->row(source));
    return seed.row(target).dot(image.row(source));
  };
  for (int string = 0; string < dimension(); ++string) {
    const double value = contract(string, string, false);
    const double derivative = directional ? contract(string, string, true) : 0.0;
    for (int orbital = 0; orbital < n_orbitals_; ++orbital)
      if ((masks_[string] >> orbital) & 1) {
        g(orbital, orbital) += value;
        dg(orbital, orbital) += derivative;
      }
  }
  // The shear topology contains every upper one-body excitation exactly once.
  // No lower generator can contribute to an upper-triangular chart derivative.
  int shear = 0;
  for (int removed = 1; removed < n_orbitals_; ++removed)
    for (int inserted = removed - 1; inserted >= 0; --inserted, ++shear)
      for (const auto& pair : shears_[shear].pairs) {
        g(inserted, removed) += pair.sign * contract(pair.target, pair.source, false);
        if (directional)
          dg(inserted, removed) += pair.sign * contract(pair.target, pair.source, true);
      }
  const Eigen::MatrixXd inverse = orbital_transform_.triangularView<Eigen::Upper>()
      .solve(Eigen::MatrixXd::Identity(n_orbitals_, n_orbitals_)).transpose();
  Eigen::MatrixXd result;
  if (directional) {
    const Eigen::MatrixXd di = -inverse * delta_upper->transpose() * inverse;
    result = dg * inverse + g * di;
  } else {
    result = g * inverse;
  }
  result.triangularView<Eigen::StrictlyLower>().setZero();
  return result;
}

std::size_t ExteriorOrbitalTransform::dynamic_bytes() const noexcept {
  std::size_t bytes =
      static_cast<std::size_t>(orbital_transform_.size()) * sizeof(double) +
      masks_.capacity() * sizeof(std::uint64_t) +
      static_cast<std::size_t>(determinant_scales_.size()) * sizeof(double) +
      shears_.capacity() * sizeof(Shear);
  for (const Shear& shear : shears_) {
    bytes += shear.pairs.capacity() * sizeof(DeterminantPair);
  }
  if (row_cache_)
    bytes += row_cache_->cached_bytes.load(std::memory_order_relaxed);
  return bytes;
}

}  // namespace xmvb::vb
