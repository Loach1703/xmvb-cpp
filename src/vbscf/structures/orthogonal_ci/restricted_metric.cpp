#include "vbscf/structures/orthogonal_ci/restricted_metric.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>

#include <Eigen/SVD>

#include "vbscf/structures/assembly/action.hpp"
#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"

namespace xmvb::vb {
namespace {

using Mask = std::uint64_t;
using Pattern = std::pair<Mask, Mask>;  // occupied union, doubly occupied

Mask occupation_mask(const std::vector<int>& string, int electrons) {
  if (static_cast<int>(string.size()) != electrons)
    throw std::invalid_argument("metric completion has inconsistent spin counts");
  Mask mask = 0;
  for (const int orbital : string) {
    if (orbital < 0 || orbital >= 63 || (mask & (Mask{1} << orbital)))
      throw std::invalid_argument("metric completion has invalid occupations");
    mask |= Mask{1} << orbital;
  }
  return mask;
}

int choose(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  std::uint64_t result = 1;
  for (int i = 1; i <= k; ++i) {
    result = result * static_cast<std::uint64_t>(n - k + i) / i;
    if (result > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
      throw std::length_error("metric spin dimension exceeds int range");
  }
  return static_cast<int>(result);
}

struct Group {
  std::vector<int> products;
  std::vector<int> selected;
  int pure_dimension = 0;
  int complement = 0;
};

bool charge_bytes(std::size_t count, std::size_t element_bytes,
                  std::size_t limit, std::size_t* total) {
  if (element_bytes == 0 || count > (limit - *total) / element_bytes)
    return false;
  *total += count * element_bytes;
  return true;
}

}  // namespace

std::unique_ptr<RestrictedMetricInverse> RestrictedMetricInverse::make(
    const StructureAction& action,
    const std::vector<std::vector<int>>& alpha_strings,
    const std::vector<std::vector<int>>& beta_strings,
    std::size_t byte_limit) {
  if (!action.storage().orthogonal_direct_ci || byte_limit == 0) return nullptr;
  try {
    return std::unique_ptr<RestrictedMetricInverse>(
        new RestrictedMetricInverse(action, alpha_strings, beta_strings, byte_limit));
  } catch (const std::invalid_argument&) {
    return nullptr;
  } catch (const std::length_error&) {
    return nullptr;
  }
}

RestrictedMetricInverse::RestrictedMetricInverse(
    const StructureAction& action,
    const std::vector<std::vector<int>>& alpha_strings,
    const std::vector<std::vector<int>>& beta_strings,
    std::size_t byte_limit)
    : n_selected_(action.n_structures_),
      n_alpha_(static_cast<int>(alpha_strings.size())),
      n_beta_(static_cast<int>(beta_strings.size())),
      byte_limit_(byte_limit) {
  if (n_alpha_ != action.n_unique_alpha_ || n_beta_ != action.n_unique_beta_ ||
      n_alpha_ <= 0 || n_beta_ <= 0 || n_selected_ <= 0 ||
      static_cast<std::size_t>(n_alpha_) * n_beta_ >
          static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("metric completion carrier dimensions disagree");
  const int n_alpha_electrons = static_cast<int>(alpha_strings.front().size());
  const int n_beta_electrons = static_cast<int>(beta_strings.front().size());
  if (n_alpha_electrons < n_beta_electrons)
    throw std::invalid_argument("metric completion requires highest-weight spin projection");

  // Preflight the full-carrier topology before allocating one entry per spin
  // product. This is a conservative bound for maps, vectors and index tables.
  const std::size_t n_products = static_cast<std::size_t>(n_alpha_) * n_beta_;
  std::size_t estimated = 0;
  if (!charge_bytes(n_products, 192, byte_limit, &estimated))
    throw std::length_error("metric completion product topology exceeds budget");
  for (const auto& s : alpha_strings)
    if (!charge_bytes(s.size(), sizeof(int), byte_limit, &estimated) ||
        !charge_bytes(1, sizeof(std::vector<int>), byte_limit, &estimated))
      throw std::length_error("metric alpha strings exceed budget");
  for (const auto& s : beta_strings)
    if (!charge_bytes(s.size(), sizeof(int), byte_limit, &estimated) ||
        !charge_bytes(1, sizeof(std::vector<int>), byte_limit, &estimated))
      throw std::length_error("metric beta strings exceed budget");
  // The inverse exterior transforms have one pair list per orbital shear.
  // Every shear contains at most one pair per spin determinant.
  for (const auto& s : alpha_strings)
    for (const int orbital : s) {
      if (orbital < 0 || orbital >= 63)
        throw std::invalid_argument("metric completion orbital index is invalid");
      n_orbitals_ = std::max(n_orbitals_, orbital + 1);
    }
  for (const auto& s : beta_strings)
    for (const int orbital : s) {
      if (orbital < 0 || orbital >= 63)
        throw std::invalid_argument("metric completion orbital index is invalid");
      n_orbitals_ = std::max(n_orbitals_, orbital + 1);
    }
  if (n_orbitals_ <= 0 ||
      choose(n_orbitals_, n_alpha_electrons) != n_alpha_ ||
      choose(n_orbitals_, n_beta_electrons) != n_beta_)
    throw std::invalid_argument("metric completion lacks a full spin-string carrier");
  const std::size_t shears = static_cast<std::size_t>(n_orbitals_) *
                             std::max(0, n_orbitals_ - 1) / 2;
  const std::size_t transform_strings =
      alpha_strings == beta_strings ? n_alpha_ : n_alpha_ + n_beta_;
  if (!charge_bytes(shears * transform_strings, 24, byte_limit, &estimated) ||
      !charge_bytes(shears, 64, byte_limit, &estimated) ||
      !charge_bytes(transform_strings, 24, byte_limit, &estimated) ||
      !charge_bytes(static_cast<std::size_t>(n_orbitals_) * n_orbitals_,
                    alpha_strings == beta_strings ? 8 : 16,
                    byte_limit, &estimated))
    throw std::length_error("metric exterior transforms exceed budget");

  alpha_strings_ = alpha_strings;
  beta_strings_ = beta_strings;

  std::vector<Mask> alpha_masks, beta_masks;
  alpha_masks.reserve(n_alpha_);
  beta_masks.reserve(n_beta_);
  for (const auto& s : alpha_strings)
    alpha_masks.push_back(occupation_mask(s, n_alpha_electrons));
  for (const auto& s : beta_strings)
    beta_masks.push_back(occupation_mask(s, n_beta_electrons));
  {
    auto unique_alpha = alpha_masks, unique_beta = beta_masks;
    std::sort(unique_alpha.begin(), unique_alpha.end());
    std::sort(unique_beta.begin(), unique_beta.end());
    if (std::adjacent_find(unique_alpha.begin(), unique_alpha.end()) != unique_alpha.end() ||
        std::adjacent_find(unique_beta.begin(), unique_beta.end()) != unique_beta.end())
      throw std::invalid_argument("metric completion has duplicate spin strings");
  }

  std::map<Pattern, Group> groups;
  std::vector<Pattern> product_patterns(static_cast<std::size_t>(n_alpha_) * n_beta_);
  for (int a = 0; a < n_alpha_; ++a)
    for (int b = 0; b < n_beta_; ++b) {
      const int product = a * n_beta_ + b;
      const Pattern pattern{alpha_masks[a] | beta_masks[b],
                            alpha_masks[a] & beta_masks[b]};
      product_patterns[product] = pattern;
      groups[pattern].products.push_back(product);
    }

  for (int s = 0; s < n_selected_; ++s) {
    const auto first = action.structure_term_offsets_[s];
    const auto last = action.structure_term_offsets_[s + 1];
    if (first == last)
      throw std::invalid_argument("metric completion has an empty selected column");
    const Pattern pattern = product_patterns.at(
        action.spin_terms_[action.structure_term_indices_[first]].spin_product);
    for (auto i = first + 1; i < last; ++i)
      if (product_patterns.at(
              action.spin_terms_[action.structure_term_indices_[i]].spin_product) != pattern)
        throw std::invalid_argument("selected structure crosses occupation blocks");
    groups.at(pattern).selected.push_back(s);
  }

  // Admission precedes any dense nullspace or complement-column allocation.
  // The bound includes the sparse incidence, block factors, transforms and a
  // largest-block decomposition workspace. One carrier RHS must also fit.
  std::size_t q_nonzeros_bound = 0;
  int max_product_block = 0;
  for (auto& [pattern, group] : groups) {
    const int n_open = __builtin_popcountll(pattern.first ^ pattern.second);
    const int n_doubly = __builtin_popcountll(pattern.second);
    const int n_alpha_open = n_alpha_electrons - n_doubly;
    const int determinant_dimension = choose(n_open, n_alpha_open);
    const int raised_dimension = choose(n_open, n_alpha_open + 1);
    if (determinant_dimension != static_cast<int>(group.products.size()))
      throw std::invalid_argument("metric completion carrier omits spin products");
    group.pure_dimension = determinant_dimension - raised_dimension;
    if (group.pure_dimension < static_cast<int>(group.selected.size()))
      throw std::invalid_argument("selected structures exceed their spin sector");
    group.complement = group.pure_dimension - static_cast<int>(group.selected.size());
    if (group.pure_dimension > std::numeric_limits<int>::max() - n_full_)
      throw std::length_error("completed spin basis exceeds int range");
    n_full_ += group.pure_dimension;
    max_product_block = std::max(max_product_block, determinant_dimension);
    const std::size_t d = group.products.size();
    const std::size_t k = static_cast<std::size_t>(group.pure_dimension);
    const std::size_t q = static_cast<std::size_t>(group.complement);
    if (q > 0 && d > (byte_limit - estimated) /
                          (sizeof(Eigen::Triplet<double, int>) + sizeof(double) + sizeof(int)) / q)
      throw std::length_error("metric complement incidence exceeds budget");
    q_nonzeros_bound += d * q;
    if (!charge_bytes(k * k, sizeof(double), byte_limit, &estimated))
      throw std::length_error("metric spin factors exceed budget");
    if (!charge_bytes(k, sizeof(int), byte_limit, &estimated) ||
        !charge_bytes(1, sizeof(Block), byte_limit, &estimated))
      throw std::length_error("metric block metadata exceeds budget");
  }
  const int n_q = n_full_ - n_selected_;
  if (n_q <= 0 || n_q >= n_selected_)
    throw std::invalid_argument("metric complement is not economical");
  const std::size_t n_nonzero_bound = action.spin_terms_.size() + q_nonzeros_bound;
  if (!charge_bytes(n_nonzero_bound,
                    sizeof(Eigen::Triplet<double, int>) + sizeof(double) + sizeof(int),
                    byte_limit, &estimated) ||
      !charge_bytes(static_cast<std::size_t>(n_full_) + 1, sizeof(int), byte_limit, &estimated) ||
      !charge_bytes(static_cast<std::size_t>(max_product_block) * max_product_block,
                    16 * sizeof(double), byte_limit, &estimated))
    throw std::length_error("metric completion exceeds its byte budget");
  const std::size_t one_column =
      (3ULL * n_products + 8ULL * n_full_) * sizeof(double);
  chunk_width_ = static_cast<int>(std::min<std::size_t>(
      8, (byte_limit - estimated) / one_column));
  if (chunk_width_ < 1)
    throw std::length_error("metric completion cannot hold one carrier RHS");

  std::vector<int> product_local(product_patterns.size(), -1);
  std::vector<Eigen::Triplet<double, int>> triplets;
  triplets.reserve(action.spin_terms_.size() + q_nonzeros_bound);
  for (const auto& term : action.spin_terms_)
    triplets.emplace_back(term.spin_product, term.structure, term.coefficient);

  int next_q = n_selected_;
  for (auto& [pattern, group] : groups) {
    const int d = static_cast<int>(group.products.size());
    const int p = static_cast<int>(group.selected.size());
    const int q = group.complement;
    const int raised_dimension = d - group.pure_dimension;
    for (int row = 0; row < d; ++row)
      product_local[group.products[row]] = row;
    Eigen::MatrixXd selected = Eigen::MatrixXd::Zero(d, p);
    for (int col = 0; col < p; ++col) {
      const int structure = group.selected[col];
      for (auto i = action.structure_term_offsets_[structure];
           i < action.structure_term_offsets_[structure + 1]; ++i) {
        const auto& term = action.spin_terms_[action.structure_term_indices_[i]];
        selected(product_local[term.spin_product], col) += term.coefficient;
      }
    }

    std::map<std::pair<Mask, Mask>, int> raised_rows;
    std::vector<std::tuple<int, int, double>> raising_terms;
    for (int col = 0; col < d; ++col) {
      const int product = group.products[col];
      const int a = product / n_beta_, b = product % n_beta_;
      for (int beta_pos = 0; beta_pos < n_beta_electrons; ++beta_pos) {
        const int orbital = beta_strings[b][beta_pos];
        const Mask bit = Mask{1} << orbital;
        if (alpha_masks[a] & bit) continue;
        const Mask lower = bit - 1;
        const int alpha_pos = __builtin_popcountll(alpha_masks[a] & lower);
        const auto key = std::make_pair(alpha_masks[a] | bit, beta_masks[b] ^ bit);
        const auto [it, inserted] = raised_rows.emplace(key,
            static_cast<int>(raised_rows.size()));
        static_cast<void>(inserted);
        const double sign = (n_alpha_electrons + beta_pos + alpha_pos) % 2 ? -1.0 : 1.0;
        raising_terms.emplace_back(it->second, col, sign);
      }
    }
    if (static_cast<int>(raised_rows.size()) != raised_dimension)
      throw std::invalid_argument("metric spin-raising sector is incomplete");
    Eigen::MatrixXd raising = Eigen::MatrixXd::Zero(raised_dimension, d);
    for (const auto& [row, col, value] : raising_terms)
      raising(row, col) += value;
    const double tolerance = 256.0 * std::numeric_limits<double>::epsilon() *
                             std::max(1, d);
    if (raised_dimension > 0) {
      const Eigen::JacobiSVD<Eigen::MatrixXd> check(raising);
      if (check.singularValues()[raised_dimension - 1] <=
          tolerance * check.singularValues()[0])
        throw std::invalid_argument("spin-raising constraints lack expected rank");
      if (p > 0 && (raising * selected).norm() >
                       tolerance * std::max(1.0, selected.norm()) *
                           std::max(1.0, raising.norm()))
        throw std::invalid_argument("selected structure is not highest-weight spin pure");
    }

    Eigen::MatrixXd normalized_selected = selected;
    for (int col = 0; col < p; ++col) {
      const double norm = selected.col(col).norm();
      if (!(norm > 0.0) || !std::isfinite(norm))
        throw std::invalid_argument("selected spin column has zero or invalid norm");
      normalized_selected.col(col) /= norm;
    }
    Eigen::MatrixXd constraints(raised_dimension + p, d);
    if (raised_dimension > 0) constraints.topRows(raised_dimension) = raising;
    if (p > 0) constraints.bottomRows(p) = normalized_selected.transpose();
    Eigen::MatrixXd missing(d, q);
    if (q > 0) {
      if (constraints.rows() == 0) {
        missing.setIdentity();
      } else {
        const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
            constraints, Eigen::ComputeFullV);
        if (svd.singularValues()[constraints.rows() - 1] <=
            tolerance * svd.singularValues()[0])
          throw std::invalid_argument("selected spin columns are linearly dependent");
        missing = svd.matrixV().rightCols(q);
      }
    } else if (p > 0) {
      const Eigen::JacobiSVD<Eigen::MatrixXd> svd(normalized_selected);
      if (svd.singularValues()[p - 1] <=
          tolerance * svd.singularValues()[0])
        throw std::invalid_argument("selected spin columns are linearly dependent");
    }

    Eigen::MatrixXd full(d, group.pure_dimension);
    if (p > 0) full.leftCols(p) = selected;
    if (q > 0) full.rightCols(q) = missing;
    const Eigen::MatrixXd gram = full.transpose() * full;
    Eigen::LLT<Eigen::MatrixXd> factor(gram);
    if (factor.info() != Eigen::Success || !gram.allFinite())
      throw std::runtime_error("completed spin metric is not positive definite");
    Block block;
    block.indices = group.selected;
    block.indices.reserve(group.pure_dimension);
    for (int col = 0; col < q; ++col) {
      const int structure = next_q++;
      block.indices.push_back(structure);
      for (int row = 0; row < d; ++row)
        if (missing(row, col) != 0.0)
          triplets.emplace_back(group.products[row], structure, missing(row, col));
    }
    block.lower = factor.matrixL();
    max_block_size_ = std::max(max_block_size_, group.pure_dimension);
    blocks_.push_back(std::move(block));
  }
  if (next_q != n_full_)
    throw std::logic_error("metric completion column count disagrees");
  expansion_.resize(n_alpha_ * n_beta_, n_full_);
  expansion_.setFromTriplets(triplets.begin(), triplets.end());
  expansion_.makeCompressed();
  stored_bytes_ = static_cast<std::size_t>(expansion_.nonZeros()) *
                      (sizeof(double) + sizeof(int)) +
                  static_cast<std::size_t>(n_full_ + 1) * sizeof(int);
  for (const auto& s : alpha_strings_) stored_bytes_ += s.size() * sizeof(int) + sizeof(s);
  for (const auto& s : beta_strings_) stored_bytes_ += s.size() * sizeof(int) + sizeof(s);
  for (const auto& block : blocks_)
    stored_bytes_ += static_cast<std::size_t>(block.lower.size()) * sizeof(double) +
                     block.indices.size() * sizeof(int);
}

Eigen::MatrixXd RestrictedMetricInverse::solve_spin(
    const Eigen::Ref<const Eigen::MatrixXd>& rhs) const {
  Eigen::MatrixXd result(rhs.rows(), rhs.cols());
  Eigen::MatrixXd local(max_block_size_, rhs.cols());
  for (const auto& block : blocks_) {
    const int d = static_cast<int>(block.indices.size());
    auto values = local.topRows(d);
    for (int row = 0; row < d; ++row) values.row(row) = rhs.row(block.indices[row]);
    block.lower.triangularView<Eigen::Lower>().solveInPlace(values);
    block.lower.transpose().triangularView<Eigen::Upper>().solveInPlace(values);
    for (int row = 0; row < d; ++row) result.row(block.indices[row]) = values.row(row);
  }
  return result;
}

Eigen::MatrixXd RestrictedMetricInverse::apply_full(
    const Eigen::Ref<const Eigen::MatrixXd>& rhs,
    const ExteriorOrbitalTransform& alpha,
    const ExteriorOrbitalTransform& beta) const {
  if (rhs.rows() != n_full_ || rhs.cols() <= 0 ||
      rhs.cols() > chunk_width_ || !rhs.allFinite())
    throw std::invalid_argument("completed metric inverse RHS is invalid");
  const Eigen::MatrixXd spin = solve_spin(rhs);
  const Eigen::MatrixXd expanded = expansion_ * spin;
  Eigen::MatrixXd carrier(n_alpha_, n_beta_ * rhs.cols());
  for (int col = 0; col < rhs.cols(); ++col)
    for (int a = 0; a < n_alpha_; ++a)
      for (int b = 0; b < n_beta_; ++b)
        carrier(a, col * n_beta_ + b) = expanded(a * n_beta_ + b, col);
  alpha.apply_inverse_adjoint_left(&carrier);
  for (int col = 0; col < rhs.cols(); ++col)
    beta.apply_inverse_adjoint_right_block(
        carrier.middleCols(col * n_beta_, n_beta_));
  alpha.apply_inverse_left(&carrier);
  for (int col = 0; col < rhs.cols(); ++col)
    beta.apply_inverse_right_block(carrier.middleCols(col * n_beta_, n_beta_));
  Eigen::MatrixXd projected(n_alpha_ * n_beta_, rhs.cols());
  for (int col = 0; col < rhs.cols(); ++col)
    for (int a = 0; a < n_alpha_; ++a)
      for (int b = 0; b < n_beta_; ++b)
        projected(a * n_beta_ + b, col) = carrier(a, col * n_beta_ + b);
  return solve_spin(expansion_.transpose() * projected);
}

RestrictedMetricInverse::Prepared::Prepared(const RestrictedMetricInverse* owner)
    : owner_(owner) {}

std::shared_ptr<const RestrictedMetricInverse::Prepared>
RestrictedMetricInverse::prepare(
    const Eigen::Ref<const Eigen::MatrixXd>& upper) const {
  if (upper.rows() != n_orbitals_ || upper.cols() != n_orbitals_)
    throw std::invalid_argument("metric trial chart differs from the complete carrier");
  auto prepared = std::shared_ptr<Prepared>(new Prepared(this));
  prepared->alpha_ = std::make_shared<ExteriorOrbitalTransform>(
      alpha_strings_, upper);
  if (alpha_strings_ == beta_strings_)
    prepared->beta_ = prepared->alpha_;
  else
    prepared->beta_ = std::make_shared<ExteriorOrbitalTransform>(
        beta_strings_, upper);
  std::size_t live_bytes = 0;
  const std::size_t transforms = prepared->alpha_->dynamic_bytes() +
      (prepared->alpha_ == prepared->beta_ ? 0 :
       prepared->beta_->dynamic_bytes());
  if (!charge_bytes(stored_bytes_, 1, byte_limit_, &live_bytes) ||
      !charge_bytes(transforms, 1, byte_limit_, &live_bytes) ||
      !charge_bytes(static_cast<std::size_t>(chunk_width_) *
                        (3ULL * n_alpha_ * n_beta_ + 8ULL * n_full_),
                    sizeof(double), byte_limit_, &live_bytes))
    throw std::length_error("prepared restricted metric exceeds its byte budget");
  return prepared;
}

Eigen::MatrixXd RestrictedMetricInverse::Prepared::precondition(
    const Eigen::Ref<const Eigen::MatrixXd>& rhs) const {
  if (rhs.rows() != owner_->n_selected_ || rhs.cols() <= 0 || !rhs.allFinite())
    throw std::invalid_argument("restricted metric RHS is invalid");
  Eigen::MatrixXd result(rhs.rows(), rhs.cols());
  for (int first = 0; first < rhs.cols(); first += owner_->chunk_width_) {
    const int count = std::min<int>(owner_->chunk_width_, rhs.cols() - first);
    Eigen::MatrixXd extended = Eigen::MatrixXd::Zero(owner_->n_full_, count);
    extended.topRows(owner_->n_selected_) = rhs.middleCols(first, count);
    const Eigen::MatrixXd image = owner_->apply_full(
        extended, *alpha_, *beta_);
    result.middleCols(first, count) = image.topRows(owner_->n_selected_);
  }
  if (!result.allFinite())
    throw std::runtime_error("restricted metric preconditioner produced nonfinite coefficients");
  return result;
}

int RestrictedMetricInverse::Prepared::n_complement() const noexcept {
  return owner_->n_complement();
}

std::size_t RestrictedMetricInverse::Prepared::stored_bytes() const noexcept {
  return owner_->stored_bytes_ + alpha_->dynamic_bytes() +
         (alpha_ == beta_ ? 0 : beta_->dynamic_bytes());
}

}  // namespace xmvb::vb
