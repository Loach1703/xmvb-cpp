#include "vbscf/structures/assembly/action.hpp"

#include <algorithm>
#include <stdexcept>

#include "core/openmp.hpp"
#include "vbscf/structures/assembly/sources.hpp"
#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"
#include "vbscf/structures/overlap/representation.hpp"

namespace xmvb::vb {
namespace {

/** @brief A bounded single-spin factor and, when requested, its tangent. */
struct Factor {
  Eigen::MatrixXd value, tangent;
  Factor(int rows, int columns, bool directional)
      : value(Eigen::MatrixXd::Zero(rows, columns)) {
    if (directional) tangent = value;
  }
};

void transform(const ExteriorOrbitalTransform& map,
    const ExteriorTransformDirection* direction, Factor& factor, bool adjoint = false) {
  if (direction) {
    if (adjoint) map.apply_directional_adjoint_left(*direction, &factor.value, &factor.tangent);
    else map.apply_directional_left(*direction, &factor.value, &factor.tangent);
  } else {
    if (adjoint) map.apply_adjoint_left(&factor.value);
    else map.apply_left(&factor.value);
  }
}

int panel_width(std::size_t bytes, int na, int nb, int columns, int copies) {
  const std::size_t row_bytes = copies * sizeof(double) *
      (static_cast<std::size_t>(na) + nb);
  if (bytes < row_bytes) throw std::length_error("cross overlap budget cannot hold one factor column");
  return static_cast<int>(std::max<std::size_t>(1, std::min<std::size_t>(columns,
      std::min(bytes / row_bytes, static_cast<std::size_t>(effective_openmp_thread_count()) * 8))));
}

} // namespace

Eigen::MatrixXd StructureAction::cross_overlap(
    const Eigen::Ref<const Eigen::MatrixXd>& vectors,
    const ExteriorOrbitalTransform& la, const ExteriorOrbitalTransform& lb,
    const ExteriorOrbitalTransform& ra, const ExteriorOrbitalTransform& rb,
    const ExteriorTransformDirection* dla, const ExteriorTransformDirection* dlb,
    const ExteriorTransformDirection* dra, const ExteriorTransformDirection* drb,
    std::size_t workspace_bytes) const {
  const bool directional = dla || dlb || dra || drb;
  if (vectors.rows() != n_structures_ || vectors.cols() <= 0 ||
      la.dimension() != n_unique_alpha_ || ra.dimension() != n_unique_alpha_ ||
      lb.dimension() != n_unique_beta_ || rb.dimension() != n_unique_beta_ ||
      (directional && !(dla && dlb && dra && drb)))
    throw std::invalid_argument("invalid selected cross overlap");
  const auto& sources = *sources_;
  const int extent = panel_width(workspace_bytes, n_unique_alpha_, n_unique_beta_,
      sources.columns.size(), directional ? 8 : 4);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n_structures_, vectors.cols());
  for (int column = 0; column < vectors.cols(); ++column) {
    Eigen::VectorXd coefficients = Eigen::VectorXd::Zero(spin_products_.size());
    Eigen::VectorXd products = Eigen::VectorXd::Zero(spin_products_.size());
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
    for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group)
      for (std::size_t term = spin_term_offsets_[group]; term < spin_term_offsets_[group + 1]; ++term)
        coefficients[group] += spin_terms_[term].coefficient * vectors(spin_terms_[term].structure, column);

    for (int begin = 0; begin < static_cast<int>(sources.columns.size()); begin += extent) {
      const int count = std::min<int>(extent, sources.columns.size() - begin);
      Factor a(n_unique_alpha_, count, directional), b(n_unique_beta_, count, directional);
      for (int j = 0; j < count; ++j) {
        const int beta = sources.columns[begin + j];
        b.value(beta, j) = 1.0;
        for (const int group : sources.beta_groups[beta])
          a.value(spin_products_[group] / n_unique_beta_, j) = coefficients[group];
      }
      // K_alpha X_J and K_beta E_J, with K_sigma = T_left^T T_right.
      // Associativity removes the transformed alpha-by-beta product entirely.
      transform(ra, dra, a); transform(la, dla, a, true);
      transform(rb, drb, b); transform(lb, dlb, b, true);
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
      for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group) {
        const int alpha = spin_products_[group] / n_unique_beta_;
        const int beta = spin_products_[group] % n_unique_beta_;
        products[group] += directional
            ? a.tangent.row(alpha).dot(b.value.row(beta)) +
                  a.value.row(alpha).dot(b.tangent.row(beta))
            : a.value.row(alpha).dot(b.value.row(beta));
      }
    }
#pragma omp parallel for schedule(static) if(n_structures_ > 32)
    for (int structure = 0; structure < n_structures_; ++structure)
      for (std::size_t offset = structure_term_offsets_[structure];
           offset < structure_term_offsets_[structure + 1]; ++offset) {
        const auto& term = spin_terms_[structure_term_indices_[offset]];
        result(structure, column) += term.coefficient * products[term.group];
      }
  }
  return result;
}

UpperPair StructureAction::cross_overlap_adjoint(
    const Eigen::VectorXd& left, const Eigen::VectorXd& right,
    const ExteriorOrbitalTransform& la, const ExteriorOrbitalTransform& lb,
    const ExteriorOrbitalTransform& ra, const ExteriorOrbitalTransform& rb,
    const Eigen::MatrixXd* delta_left, const Eigen::MatrixXd* delta_right,
    std::size_t workspace_bytes) const {
  const bool directional = delta_left != nullptr;
  if (left.size() != n_structures_ || right.size() != n_structures_ ||
      directional != (delta_right != nullptr))
    throw std::invalid_argument("invalid cross overlap adjoint");
  const auto& sources = *sources_;
  const int extent = panel_width(workspace_bytes, n_unique_alpha_, n_unique_beta_,
      sources.columns.size(), directional ? 16 : 8);
  ExteriorTransformDirection dla, dlb, dra, drb;
  if (directional) {
    dla = la.direction(*delta_left); dlb = lb.direction(*delta_left);
    dra = ra.direction(*delta_right); drb = rb.direction(*delta_right);
  }
  Eigen::MatrixXd coefficients = Eigen::MatrixXd::Zero(spin_products_.size(), 2);
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
  for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group)
    for (std::size_t t = spin_term_offsets_[group]; t < spin_term_offsets_[group + 1]; ++t) {
      const auto& term = spin_terms_[t];
      coefficients(group, 0) += term.coefficient * left[term.structure];
      coefficients(group, 1) += term.coefficient * right[term.structure];
    }
  UpperPair result{Eigen::MatrixXd::Zero(la.n_orbitals(), la.n_orbitals()),
                   Eigen::MatrixXd::Zero(ra.n_orbitals(), ra.n_orbitals())};
  for (int begin = 0; begin < static_cast<int>(sources.columns.size()); begin += extent) {
    const int count = std::min<int>(extent, sources.columns.size() - begin);
    Factor u(n_unique_alpha_, count, directional), v(n_unique_beta_, count, directional);
    for (int j = 0; j < count; ++j) {
      const int beta = sources.columns[begin + j];
      v.value(beta, j) = 1.0;
      for (const int group : sources.beta_groups[beta])
        u.value(spin_products_[group] / n_unique_beta_, j) = coefficients(group, 1);
    }
    transform(ra, directional ? &dra : nullptr, u);
    transform(rb, directional ? &drb : nullptr, v);
    Factor a = u, b = v;
    transform(la, directional ? &dla : nullptr, a, true);
    transform(lb, directional ? &dlb : nullptr, b, true);
    Factor za(n_unique_alpha_, count, directional), zb(n_unique_beta_, count, directional);
    // Seeds X_left B and X_left^T A use only the original touched products.
#pragma omp parallel for schedule(static) if(n_unique_alpha_ > 32)
    for (int alpha = 0; alpha < n_unique_alpha_; ++alpha)
      for (int group = sources.alpha_offsets[alpha]; group < sources.alpha_offsets[alpha + 1]; ++group) {
        const int beta = spin_products_[group] % n_unique_beta_;
        za.value.row(alpha).noalias() += coefficients(group, 0) * b.value.row(beta);
        if (directional) za.tangent.row(alpha).noalias() += coefficients(group, 0) * b.tangent.row(beta);
      }
#pragma omp parallel for schedule(static) if(n_unique_beta_ > 32)
    for (int beta = 0; beta < n_unique_beta_; ++beta)
      for (const int group : sources.beta_groups[beta]) {
        const int alpha = spin_products_[group] / n_unique_beta_;
        zb.value.row(beta).noalias() += coefficients(group, 0) * a.value.row(alpha);
        if (directional) zb.tangent.row(beta).noalias() += coefficients(group, 0) * a.tangent.row(alpha);
      }
    transform(la, directional ? &dla : nullptr, za);
    transform(lb, directional ? &dlb : nullptr, zb);
    const auto pullback = [&](const ExteriorOrbitalTransform& map,
        const Factor& seed, const Factor& image, const Eigen::MatrixXd* delta) {
      return map.generator_pullback(seed.value, image.value,
          directional ? &seed.tangent : nullptr, directional ? &image.tangent : nullptr, delta);
    };
    result.left += pullback(la, u, za, delta_left) + pullback(lb, v, zb, delta_left);
    result.right += pullback(ra, za, u, delta_right) + pullback(rb, zb, v, delta_right);
  }
  return result;
}

} // namespace xmvb::vb
