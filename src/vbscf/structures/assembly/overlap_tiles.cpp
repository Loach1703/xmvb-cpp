#include "vbscf/structures/assembly/action.hpp"

#include <algorithm>
#include <stdexcept>

#include "core/openmp.hpp"
#include "mcscf/ci/product_support.hpp"
#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"
#include "vbscf/structures/overlap/representation.hpp"

namespace xmvb::vb {
namespace {

/** @brief Selected rows of d(wedge R), without a dense compound matrix. */
Eigen::MatrixXd directional_rows(const ExteriorOrbitalTransform& transform,
    const ExteriorTransformDirection& direction, int begin, int count) {
  Eigen::MatrixXd rows = Eigen::MatrixXd::Zero(count, transform.dimension());
  rows.middleCols(begin, count).setIdentity();
  Eigen::MatrixXd derivative = Eigen::MatrixXd::Zero(count, transform.dimension());
  transform.apply_directional_adjoint_right(direction, &rows, &derivative);
  return derivative;
}

/** @brief Panels intersecting exact triangular support, without a product bitmap. */
std::vector<int> reachable_panels(const mcscf::ProductSupport* support,
    int n_beta, int extent) {
  std::vector<int> panels;
  if (support) {
    auto cursor = support->cursor(0);
    int beta;
    while (cursor.next_column(beta)) panels.push_back(beta / extent);
    std::sort(panels.begin(), panels.end());
    panels.erase(std::unique(panels.begin(), panels.end()), panels.end());
  } else {
    for (int first = 0; first < n_beta; first += extent) panels.push_back(first / extent);
  }
  return panels;
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

  // Budget counts concurrent primal/tangent panels and beta compound rows.
  // At least one column must fit; no change of algorithm or physical space.
  const std::size_t copies = directional ? 8 : 4;
  const std::size_t row_bytes = copies * sizeof(double) *
      (static_cast<std::size_t>(n_unique_alpha_) + n_unique_beta_);
  if (workspace_bytes < row_bytes)
    throw std::length_error("cross overlap budget cannot hold one column");
  const int extent = static_cast<int>(std::min<std::size_t>(n_unique_beta_,
      std::min(workspace_bytes / row_bytes,
          static_cast<std::size_t>(effective_openmp_thread_count()) * 8)));

  // Generate reachable columns, rather than visiting the whole product grid.
  // Upper-triangular chart derivatives preserve this structural support.
  const auto panels = reachable_panels(product_support_.get(), n_unique_beta_, extent);

  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n_structures_, vectors.cols());
  // One RHS at a time bounds storage independently of Davidson block width.
  for (int column = 0; column < vectors.cols(); ++column) {
    Eigen::VectorXd coefficients = Eigen::VectorXd::Zero(spin_products_.size());
    Eigen::VectorXd products = Eigen::VectorXd::Zero(spin_products_.size());
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
    for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group)
      for (std::size_t term = spin_term_offsets_[group];
           term < spin_term_offsets_[group + 1]; ++term)
        coefficients[group] += spin_terms_[term].coefficient *
            vectors(spin_terms_[term].structure, column);

    for (const int panel : panels) {
      const int begin = panel * extent;
      const int count = std::min(extent, n_unique_beta_ - begin);
      const auto right_rows = rb.row_block(begin, count);
      const auto left_rows = (&lb == &rb) ? right_rows : lb.row_block(begin, count);
      Eigen::MatrixXd d_right, d_left;
      if (directional) {
        d_right = directional_rows(rb, *drb, begin, count);
        d_left = directional_rows(lb, *dlb, begin, count);
      }
      Eigen::MatrixXd image = Eigen::MatrixXd::Zero(n_unique_alpha_, count);
      Eigen::MatrixXd derivative;
      if (directional) derivative = Eigen::MatrixXd::Zero(n_unique_alpha_, count);
#pragma omp parallel for schedule(static) if(n_unique_alpha_ > 32)
      for (int alpha = 0; alpha < n_unique_alpha_; ++alpha) {
        const int first = alpha * n_unique_beta_;
        auto group = std::lower_bound(spin_products_.begin(), spin_products_.end(), first);
        for (; group != spin_products_.end() && *group < first + n_unique_beta_; ++group) {
          const auto id = group - spin_products_.begin();
          const int beta = *group - first;
          image.row(alpha).noalias() += coefficients[id] * right_rows->col(beta).transpose();
          if (directional)
            derivative.row(alpha).noalias() += coefficients[id] * d_right.col(beta).transpose();
        }
      }
      if (directional && product_support_) {
        ra.apply_supported_left(&image, *product_support_, begin, count, false, dra, &derivative);
        la.apply_supported_left(&image, *product_support_, begin, count, true, dla, &derivative);
      } else if (directional) {
        ra.apply_directional_left(*dra, &image, &derivative);
        la.apply_directional_adjoint_left(*dla, &image, &derivative);
      } else if (product_support_) {
        ra.apply_supported_left(&image, *product_support_, begin, count);
        la.apply_supported_left(&image, *product_support_, begin, count, true);
      } else {
        ra.apply_left(&image);
        la.apply_adjoint_left(&image);
      }
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
      for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group) {
        const int alpha = spin_products_[group] / n_unique_beta_;
        const int beta = spin_products_[group] % n_unique_beta_;
        products[group] += directional
            ? derivative.row(alpha).dot(left_rows->col(beta)) +
                  image.row(alpha).dot(d_left.col(beta))
            : image.row(alpha).dot(left_rows->col(beta));
      }
    }
    // One final structure projection, after all panels have been consumed.
#pragma omp parallel for schedule(static) if(n_structures_ > 32)
    for (int structure = 0; structure < n_structures_; ++structure)
      for (std::size_t offset = structure_term_offsets_[structure];
           offset < structure_term_offsets_[structure + 1]; ++offset) {
        const auto& term = spin_terms_[structure_term_indices_[offset]];
        const auto group = std::lower_bound(spin_products_.begin(), spin_products_.end(),
            term.spin_product) - spin_products_.begin();
        result(structure, column) += term.coefficient * products[group];
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
  const std::size_t row_bytes = (directional ? 16 : 8) * sizeof(double) *
      (static_cast<std::size_t>(n_unique_alpha_) + n_unique_beta_);
  if (workspace_bytes < row_bytes)
    throw std::length_error("cross overlap adjoint budget cannot hold one column");
  const int extent = static_cast<int>(std::min<std::size_t>(n_unique_beta_,
      std::min(workspace_bytes / row_bytes,
          static_cast<std::size_t>(effective_openmp_thread_count()) * 8)));
  const auto panels = reachable_panels(product_support_.get(), n_unique_beta_, extent);
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
  // This index is linear in touched products, never in the ambient grid.
  std::vector<std::vector<int>> beta_groups(n_unique_beta_);
  std::vector<std::uint8_t> beta_seed(n_unique_beta_, 0);
  for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group)
    beta_groups[spin_products_[group] % n_unique_beta_].push_back(group);
  for (int beta = 0; beta < n_unique_beta_; ++beta)
    beta_seed[beta] = !beta_groups[beta].empty();

  UpperPair result{Eigen::MatrixXd::Zero(la.n_orbitals(), la.n_orbitals()),
                   Eigen::MatrixXd::Zero(ra.n_orbitals(), ra.n_orbitals())};
  for (const int panel : panels) {
    const int begin = panel * extent;
    const int count = std::min(extent, n_unique_beta_ - begin);
    const auto lr = lb.row_block(begin, count), rr = rb.row_block(begin, count);
    Eigen::MatrixXd dlr, drr;
    if (directional) {
      dlr = directional_rows(lb, dlb, begin, count);
      drr = directional_rows(rb, drb, begin, count);
    }
    Eigen::MatrixXd yl, yr, dyl, dyr;
    const auto expand = [&](int side, const Eigen::MatrixXd& rows,
        const Eigen::MatrixXd& drows, const ExteriorOrbitalTransform& alpha,
        const ExteriorTransformDirection& da, Eigen::MatrixXd& y, Eigen::MatrixXd& dy) {
      y = Eigen::MatrixXd::Zero(n_unique_alpha_, count);
      if (directional) dy = y;
#pragma omp parallel for schedule(static) if(n_unique_alpha_ > 32)
      for (int a = 0; a < n_unique_alpha_; ++a) {
        const int first = a * n_unique_beta_;
        auto group = std::lower_bound(spin_products_.begin(), spin_products_.end(), first);
        for (; group != spin_products_.end() && *group < first + n_unique_beta_; ++group) {
          const double c = coefficients(group - spin_products_.begin(), side);
          y.row(a).noalias() += c * rows.col(*group - first).transpose();
          if (directional) dy.row(a).noalias() += c * drows.col(*group - first).transpose();
        }
      }
      if (product_support_)
        alpha.apply_supported_left(&y, *product_support_, begin, count, false,
            directional ? &da : nullptr, directional ? &dy : nullptr);
      else if (directional) alpha.apply_directional_left(da, &y, &dy);
      else alpha.apply_left(&y);
    };
    expand(0, *lr, dlr, la, dla, yl, dyl);
    expand(1, *rr, drr, ra, dra, yr, dyr);
    result.left += la.generator_pullback(yr, yl, false,
        directional ? &dyr : nullptr, directional ? &dyl : nullptr, delta_left);
    result.right += ra.generator_pullback(yl, yr, false,
        directional ? &dyl : nullptr, directional ? &dyr : nullptr, delta_right);

    const auto beta_pullback = [&](int side, const ExteriorOrbitalTransform& alpha,
        const ExteriorOrbitalTransform& beta, const ExteriorTransformDirection& da,
        Eigen::MatrixXd& y, Eigen::MatrixXd& dy, const Eigen::MatrixXd& rows,
        const Eigen::MatrixXd& drows, const Eigen::MatrixXd* delta) {
      if (product_support_)
        alpha.apply_supported_left(&y, *product_support_, begin, count, true,
            directional ? &da : nullptr, directional ? &dy : nullptr);
      else if (directional) alpha.apply_directional_adjoint_left(da, &y, &dy);
      else alpha.apply_adjoint_left(&y);
      Eigen::MatrixXd seed = Eigen::MatrixXd::Zero(count, n_unique_beta_);
      Eigen::MatrixXd dseed;
      if (directional) dseed = seed;
#pragma omp parallel for schedule(static) if(n_unique_beta_ > 32)
      for (int b = 0; b < n_unique_beta_; ++b)
        for (int group : beta_groups[b]) {
          const int a = spin_products_[group] / n_unique_beta_;
          seed.col(b).noalias() += coefficients(group, side) * y.row(a).transpose();
          if (directional)
            dseed.col(b).noalias() += coefficients(group, side) * dy.row(a).transpose();
        }
      return beta.generator_pullback(seed, rows, true,
          directional ? &dseed : nullptr, directional ? &drows : nullptr, delta, &beta_seed);
    };
    // Both alpha pullbacks have consumed the primal images. Reuse their
    // buffers for the two beta pullbacks instead of retaining another carrier.
    result.left += beta_pullback(0, la, lb, dla, yr, dyr, *lr, dlr, delta_left);
    result.right += beta_pullback(1, ra, rb, dra, yl, dyl, *rr, drr, delta_right);
  }
  return result;
}

} // namespace xmvb::vb
