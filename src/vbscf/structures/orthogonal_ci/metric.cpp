#include "vbscf/structures/orthogonal_ci/metric.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

#include <lapacke.h>

#include "vbscf/structures/assembly/action.hpp"

namespace xmvb::vb {
namespace {

std::uint64_t occupation_mask(const std::vector<int>& string) {
  std::uint64_t mask = 0;
  for (int orbital : string) {
    if (orbital < 0 || orbital >= 63 || (mask & (std::uint64_t{1} << orbital)) != 0)
      throw std::invalid_argument("structure metric has invalid spin occupations");
    mask |= std::uint64_t{1} << orbital;
  }
  return mask;
}

}  // namespace

StructureMetricInverse::StructureMetricInverse(
    const StructureAction& action,
    const std::vector<std::vector<int>>& alpha_strings,
    const std::vector<std::vector<int>>& beta_strings,
    bool complete_space,
    std::size_t factor_byte_limit,
    std::size_t workspace_byte_limit)
    : action_(&action) {
  if (!complete_space || !action.storage().orthogonal_direct_ci)
    throw std::invalid_argument("factorized metric inverse requires a complete spin sector");
  if (alpha_strings.size() != static_cast<std::size_t>(action.n_unique_alpha_) ||
      beta_strings.size() != static_cast<std::size_t>(action.n_unique_beta_))
    throw std::invalid_argument("structure metric spin dimensions disagree");
  const std::size_t n = action.n_structures_;
  const std::size_t n_products = alpha_strings.size() * beta_strings.size();
  const std::size_t column_bytes = (n_products + 4 * n) * sizeof(double);
  if (column_bytes == 0 || workspace_byte_limit < column_bytes)
    throw std::length_error("structure metric workspace cannot hold one carrier column");
  chunk_width_ = static_cast<int>(std::min<std::size_t>(
      workspace_byte_limit / column_bytes,
      std::numeric_limits<int>::max() / action.n_unique_beta_));

  using Pattern = std::pair<std::uint64_t, std::uint64_t>;
  std::vector<std::uint64_t> alpha_masks, beta_masks;
  for (const auto& string : alpha_strings) alpha_masks.push_back(occupation_mask(string));
  for (const auto& string : beta_strings) beta_masks.push_back(occupation_mask(string));
  const auto pattern = [&](int spin_product) {
    const auto alpha = alpha_masks.at(spin_product / action.n_unique_beta_);
    const auto beta = beta_masks.at(spin_product % action.n_unique_beta_);
    return Pattern{alpha | beta, alpha & beta};
  };
  std::map<Pattern, std::vector<int>> groups;
  for (int structure = 0; structure < action.n_structures_; ++structure) {
    const auto first = action.structure_term_offsets_[structure];
    const auto last = action.structure_term_offsets_[structure + 1];
    if (first == last) throw std::invalid_argument("structure metric expansion has an empty column");
    const Pattern occupation = pattern(action.spin_terms_[action.structure_term_indices_[first]].spin_product);
    for (auto position = first + 1; position < last; ++position) {
      if (pattern(action.spin_terms_[action.structure_term_indices_[position]].spin_product) != occupation)
        throw std::invalid_argument("structure spans multiple spatial occupation patterns");
    }
    groups[occupation].push_back(structure);
  }
  // Admit all factor storage before allocating any dense block.
  for (const auto& group : groups) {
    const std::size_t d = group.second.size();
    if (d > (factor_byte_limit - factor_bytes_) / sizeof(double) / d)
      throw std::length_error("spin-recoupling factors exceed their byte budget");
    factor_bytes_ += d * d * sizeof(double);
  }
  blocks_.reserve(groups.size());
  std::vector<int> group_id(n), local_id(n);
  for (auto& group : groups) {
    const int index = static_cast<int>(blocks_.size());
    const int d = static_cast<int>(group.second.size());
    for (int local = 0; local < d; ++local) {
      group_id[group.second[local]] = index;
      local_id[group.second[local]] = local;
    }
    blocks_.push_back({std::move(group.second), Eigen::MatrixXd::Zero(d, d)});
    max_block_size_ = std::max(max_block_size_, d);
  }
  // C_0^T C_0 from sparse spin-product rows; never form the dense C_0 or Gram.
  for (std::size_t product = 0; product < action.spin_products_.size(); ++product) {
    const auto first = action.spin_term_offsets_[product];
    const auto last = action.spin_term_offsets_[product + 1];
    for (auto left = first; left < last; ++left) {
      const auto& a = action.spin_terms_[left];
      for (auto right = first; right < last; ++right) {
        const auto& b = action.spin_terms_[right];
        if (group_id[a.structure] != group_id[b.structure])
          throw std::logic_error("one spin product touches inconsistent occupation blocks");
        const int row = local_id[a.structure];
        const int column = local_id[b.structure];
        if (row >= column)
          blocks_[group_id[a.structure]].lower(row, column) += a.coefficient * b.coefficient;
      }
    }
  }
  for (auto& block : blocks_) {
    const lapack_int d = static_cast<lapack_int>(block.lower.rows());
    const lapack_int info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', d, block.lower.data(), d);
    if (info != 0 || !block.lower.allFinite())
      throw std::runtime_error("spin-recoupling metric block is not numerically positive definite");
  }
}

Eigen::MatrixXd StructureMetricInverse::solve_spin(
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

Eigen::MatrixXd StructureMetricInverse::apply(
    const Eigen::Ref<const Eigen::MatrixXd>& rhs) const {
  if (rhs.rows() != action_->n_structures_ || rhs.cols() <= 0 || !rhs.allFinite())
    throw std::invalid_argument("structure metric inverse RHS is invalid");
  Eigen::MatrixXd result(rhs.rows(), rhs.cols());
  for (Eigen::Index first = 0; first < rhs.cols(); first += chunk_width_) {
    const Eigen::Index count = std::min<Eigen::Index>(chunk_width_, rhs.cols() - first);
    Eigen::MatrixXd carrier = action_->expand_structure_block(solve_spin(rhs.middleCols(first, count)));
    action_->inverse_adjoint_orbital_transform(&carrier);
    action_->inverse_orbital_transform(&carrier);
    result.middleCols(first, count) = solve_spin(action_->contract_spin_product_block(carrier));
  }
  if (!result.allFinite())
    throw std::runtime_error("structure metric inverse produced nonfinite coefficients");
  return result;
}

}  // namespace xmvb::vb
