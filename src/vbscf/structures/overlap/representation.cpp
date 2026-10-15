#include "vbscf/structures/overlap/representation.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

#include <Eigen/LU>

#include "vbscf/structures/assembly/action.hpp"
#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"

namespace xmvb::vb {
namespace {

using Strings = std::vector<std::vector<int>>;
constexpr std::size_t kRowWorkspaceBytes = 8 * 1024 * 1024;

void validate_strings(const Strings& strings, int orbitals) {
  if (strings.empty()) throw std::invalid_argument("empty spin-string carrier");
  const std::size_t electrons = strings.front().size();
  for (const auto& occupied : strings) {
    if (occupied.size() != electrons ||
        !std::is_sorted(occupied.begin(), occupied.end()) ||
        std::adjacent_find(occupied.begin(), occupied.end()) != occupied.end() ||
        (!occupied.empty() && (occupied.front() < 0 || occupied.back() >= orbitals))) {
      throw std::invalid_argument("invalid spin-string occupation list");
    }
  }
  std::vector<const std::vector<int>*> ordered;
  ordered.reserve(strings.size());
  for (const auto& occupied : strings) ordered.push_back(&occupied);
  std::sort(ordered.begin(), ordered.end(), [](const auto* lhs, const auto* rhs) {
    return *lhs < *rhs;
  });
  if (std::adjacent_find(ordered.begin(), ordered.end(),
                         [](const auto* lhs, const auto* rhs) {
                           return *lhs == *rhs;
                         }) != ordered.end())
    throw std::invalid_argument("duplicate spin string");
}

void validate_chart(const Eigen::Ref<const Eigen::MatrixXd>& r, int n) {
  if (r.rows() != n || r.cols() != n || !r.allFinite() ||
      (r.diagonal().array() <= 0.0).any())
    throw std::invalid_argument("invalid orbital chart");
  if (n > 1) {
    const double scale = std::max(1.0, r.cwiseAbs().maxCoeff());
    if (r.template triangularView<Eigen::StrictlyLower>()
            .toDenseMatrix().cwiseAbs().maxCoeff() >
        64.0 * std::numeric_limits<double>::epsilon() * scale)
      throw std::invalid_argument("orbital chart must be upper triangular");
  }
}

void validate_direction(const Eigen::Ref<const Eigen::MatrixXd>& d, int n) {
  if (d.rows() != n || d.cols() != n || !d.allFinite())
    throw std::invalid_argument("invalid orbital chart direction");
  if (n > 1) {
    const double scale = std::max(1.0, d.cwiseAbs().maxCoeff());
    if (d.template triangularView<Eigen::StrictlyLower>()
            .toDenseMatrix().cwiseAbs().maxCoeff() >
        64.0 * std::numeric_limits<double>::epsilon() * scale)
      throw std::invalid_argument("chart direction must be upper triangular");
  }
}

Eigen::MatrixXd upper_part(const Eigen::MatrixXd& matrix) {
  Eigen::MatrixXd result = matrix;
  result.template triangularView<Eigen::StrictlyLower>().setZero();
  return result;
}

// The determinant of a deleted minor is defined also at singular K. This is
// essential: neither the cofactor nor its derivative uses K^{-1}.
double deleted_minor(const Eigen::MatrixXd& k, const std::vector<int>& rows,
                     const std::vector<int>& cols, int omit_row_1 = -1,
                     int omit_col_1 = -1, int omit_row_2 = -1,
                     int omit_col_2 = -1) {
  const int size = static_cast<int>(rows.size()) - (omit_row_1 >= 0) -
                   (omit_row_2 >= 0);
  if (size == 0) return 1.0;
  Eigen::MatrixXd minor(size, size);
  for (int i = 0, ii = 0; i < static_cast<int>(rows.size()); ++i) {
    if (i == omit_row_1 || i == omit_row_2) continue;
    for (int j = 0, jj = 0; j < static_cast<int>(cols.size()); ++j) {
      if (j == omit_col_1 || j == omit_col_2) continue;
      minor(ii, jj++) = k(rows[i], cols[j]);
    }
    ++ii;
  }
  return minor.partialPivLu().determinant();
}

double pair_value(const Eigen::MatrixXd& k, const std::vector<int>& left,
                  const std::vector<int>& right) {
  return deleted_minor(k, left, right);
}

double pair_direction(const Eigen::MatrixXd& k, const Eigen::MatrixXd& dk,
                      const std::vector<int>& left,
                      const std::vector<int>& right) {
  double result = 0.0;
  const int electrons = static_cast<int>(left.size());
  for (int i = 0; i < electrons; ++i)
    for (int j = 0; j < electrons; ++j) {
      const double sign = (i + j) % 2 ? -1.0 : 1.0;
      result += sign * deleted_minor(k, left, right, i, j) *
                dk(left[i], right[j]);
    }
  return result;
}

void add_pair_gradient(const Eigen::MatrixXd& k, const Eigen::MatrixXd* dk,
                       const std::vector<int>& left,
                       const std::vector<int>& right, double weight,
                       double dweight, Eigen::MatrixXd* gradient,
                       Eigen::MatrixXd* dgradient) {
  const int electrons = static_cast<int>(left.size());
  for (int i = 0; i < electrons; ++i)
    for (int j = 0; j < electrons; ++j) {
      const double sign = (i + j) % 2 ? -1.0 : 1.0;
      const double cofactor = sign * deleted_minor(k, left, right, i, j);
      (*gradient)(left[i], right[j]) += weight * cofactor;
      if (dgradient == nullptr) continue;
      double dcofactor = 0.0;
      for (int p = 0; p < electrons; ++p) {
        if (p == i) continue;
        for (int q = 0; q < electrons; ++q) {
          if (q == j) continue;
          const int reduced_p = p - (p > i);
          const int reduced_q = q - (q > j);
          const double second_sign = (reduced_p + reduced_q) % 2 ? -sign : sign;
          dcofactor += second_sign *
              deleted_minor(k, left, right, i, j, p, q) *
              (*dk)(left[p], right[q]);
        }
      }
      (*dgradient)(left[i], right[j]) +=
          dweight * cofactor + weight * dcofactor;
    }
}

Eigen::MatrixXd row_of(const StructureAction& action,
                       const Eigen::MatrixXd& coeff, int alpha, int beta_size) {
  return action.expand_structure_tile(coeff, alpha, 1, 0, beta_size);
}

Eigen::MatrixXd column_of(const StructureAction& action,
                          const Eigen::MatrixXd& coeff, int beta, int alpha_size) {
  return action.expand_structure_tile(coeff, 0, alpha_size, beta, 1);
}

Eigen::MatrixXd carrier_image(const StructureAction& action,
                              const ExteriorOrbitalTransform& alpha,
                              const ExteriorOrbitalTransform& beta,
                              const Eigen::MatrixXd& coeff,
                              const ExteriorTransformDirection* da = nullptr,
                              const ExteriorTransformDirection* db = nullptr,
                              Eigen::MatrixXd* derivative = nullptr) {
  Eigen::MatrixXd image = action.expand_structure_block(coeff);
  if (derivative != nullptr) {
    *derivative = Eigen::MatrixXd::Zero(image.rows(), image.cols());
    alpha.apply_directional_left(*da, &image, derivative);
    beta.apply_directional_right(*db, &image, derivative);
  } else {
    alpha.apply_left(&image);
    beta.apply_right(&image);
  }
  return image;
}

Eigen::MatrixXd carrier_adjoint(
    const StructureAction& action, const ExteriorOrbitalTransform& alpha,
    const ExteriorOrbitalTransform& beta, Eigen::MatrixXd image,
    const ExteriorTransformDirection* da = nullptr,
    const ExteriorTransformDirection* db = nullptr,
    Eigen::MatrixXd* derivative = nullptr) {
  if (derivative != nullptr) {
    alpha.apply_directional_adjoint_left(*da, &image, derivative);
    beta.apply_directional_adjoint_right(*db, &image, derivative);
  } else {
    alpha.apply_adjoint_left(&image);
    beta.apply_adjoint_right(&image);
  }
  return action.contract_spin_product_block(image);
}

Eigen::MatrixXd apply_pair_kernel(
    const StructureAction& action, const Strings& alpha, const Strings& beta,
    const Eigen::MatrixXd& kernel,
    const Eigen::Ref<const Eigen::MatrixXd>& coeff) {
  const int na = static_cast<int>(alpha.size());
  const int nb = static_cast<int>(beta.size());
  Eigen::MatrixXd result(action.n_structures(), coeff.cols());
  const std::size_t denominator = sizeof(double) *
      static_cast<std::size_t>(std::max(1, 2 * nb + 1));
  const int rows_per_tile = std::max(1, std::min(
      na, static_cast<int>(kRowWorkspaceBytes / denominator)));
  for (int col = 0; col < coeff.cols(); ++col) {
    Eigen::MatrixXd one_coeff = coeff.col(col);
    Eigen::MatrixXd one_result = Eigen::MatrixXd::Zero(action.n_structures(), 1);
    for (int begin = 0; begin < na; begin += rows_per_tile) {
      const int rows = std::min(rows_per_tile, na - begin);
      Eigen::MatrixXd projected = Eigen::MatrixXd::Zero(rows, nb);
      for (int ar = 0; ar < na; ++ar) {
        const Eigen::MatrixXd source = row_of(action, one_coeff, ar, nb);
        for (int ai = 0; ai < rows; ++ai)
          projected.row(ai).noalias() +=
              pair_value(kernel, alpha[begin + ai], alpha[ar]) * source.row(0);
      }
      Eigen::MatrixXd tile(rows, nb);
      for (int bl = 0; bl < nb; ++bl) {
        Eigen::VectorXd overlap(nb);
        for (int br = 0; br < nb; ++br)
          overlap[br] = pair_value(kernel, beta[bl], beta[br]);
        tile.col(bl).noalias() = projected * overlap;
      }
      action.add_spin_product_tile(tile, begin, 0, &one_result);
    }
    result.col(col) = one_result;
  }
  return result;
}

Eigen::MatrixXd chart_inverse_transpose(
    const Eigen::Ref<const Eigen::MatrixXd>& r) {
  return r.template triangularView<Eigen::Upper>()
      .solve(Eigen::MatrixXd::Identity(r.rows(), r.cols())).transpose();
}

}  // namespace

RepresentationOverlap::RepresentationOverlap(
    const StructureAction& action, const Strings& alpha_strings,
    const Strings& beta_strings, int n_active)
    : action_(&action), alpha_(alpha_strings), beta_(beta_strings),
      n_active_(n_active), complete_carrier_(action.supports_integral_direction()) {
  if (n_active <= 0 || (complete_carrier_ && n_active > 63))
    throw std::invalid_argument("invalid orbital count for representation overlap");
  validate_strings(alpha_, n_active);
  validate_strings(beta_, n_active);
  if (static_cast<int>(alpha_.size()) != action.n_alpha_strings() ||
      static_cast<int>(beta_.size()) != action.n_beta_strings())
    throw std::invalid_argument("spin strings do not match action carrier");
}

RepresentationOverlap::PreparedMetric::PreparedMetric(
    const StructureAction& action, const Strings& alpha, const Strings& beta,
    bool complete_carrier, const Eigen::Ref<const Eigen::MatrixXd>& upper)
    : action_(&action), alpha_(&alpha), beta_(&beta) {
  if (complete_carrier) {
    alpha_transform_ =
        std::make_shared<ExteriorOrbitalTransform>(alpha, upper);
    if (&alpha == &beta || alpha == beta) {
      beta_transform_ = alpha_transform_;
    } else {
      beta_transform_ =
          std::make_shared<ExteriorOrbitalTransform>(beta, upper);
    }
  } else {
    kernel_ = upper.transpose() * upper;
  }
}

Eigen::MatrixXd RepresentationOverlap::PreparedMetric::operator()(
    const Eigen::Ref<const Eigen::MatrixXd>& coeff) const {
  if (coeff.rows() != action_->n_structures() || coeff.cols() <= 0 ||
      !coeff.allFinite()) throw std::invalid_argument("invalid structure block");
  if (!alpha_transform_)
    return apply_pair_kernel(*action_, *alpha_, *beta_, kernel_, coeff);

  Eigen::MatrixXd result(action_->n_structures(), coeff.cols());
  for (int col = 0; col < coeff.cols(); ++col) {
    Eigen::MatrixXd one = coeff.col(col);
    Eigen::MatrixXd image = carrier_image(
        *action_, *alpha_transform_, *beta_transform_, one);
    result.col(col) = carrier_adjoint(
        *action_, *alpha_transform_, *beta_transform_, std::move(image));
  }
  return result;
}

RepresentationOverlap::PreparedMetric RepresentationOverlap::prepare_metric(
    const Eigen::Ref<const Eigen::MatrixXd>& upper) const {
  validate_chart(upper, n_active_);
  return PreparedMetric(*action_, alpha_, beta_, complete_carrier_, upper);
}

Eigen::MatrixXd RepresentationOverlap::apply(
    const Eigen::Ref<const Eigen::MatrixXd>& left_upper,
    const Eigen::Ref<const Eigen::MatrixXd>& right_upper,
    const Eigen::Ref<const Eigen::MatrixXd>& coeff) const {
  validate_chart(left_upper, n_active_);
  validate_chart(right_upper, n_active_);
  if (coeff.rows() != action_->n_structures() || coeff.cols() <= 0 ||
      !coeff.allFinite()) throw std::invalid_argument("invalid structure block");
  Eigen::MatrixXd result(action_->n_structures(), coeff.cols());
  if (complete_carrier_) {
    ExteriorOrbitalTransform la(alpha_, left_upper), lb(beta_, left_upper);
    ExteriorOrbitalTransform ra(alpha_, right_upper), rb(beta_, right_upper);
    for (int col = 0; col < coeff.cols(); ++col) {
      Eigen::MatrixXd x = coeff.col(col);
      Eigen::MatrixXd y = carrier_image(*action_, ra, rb, x);
      result.col(col) = carrier_adjoint(*action_, la, lb, std::move(y));
    }
    return result;
  }
  return apply_pair_kernel(*action_, alpha_, beta_,
                           left_upper.transpose() * right_upper, coeff);
}

Eigen::MatrixXd RepresentationOverlap::apply_direction(
    const Eigen::Ref<const Eigen::MatrixXd>& left_upper,
    const Eigen::Ref<const Eigen::MatrixXd>& right_upper,
    const Eigen::Ref<const Eigen::MatrixXd>& delta_left_upper,
    const Eigen::Ref<const Eigen::MatrixXd>& delta_right_upper,
    const Eigen::Ref<const Eigen::MatrixXd>& coeff) const {
  validate_chart(left_upper, n_active_);
  validate_chart(right_upper, n_active_);
  validate_direction(delta_left_upper, n_active_);
  validate_direction(delta_right_upper, n_active_);
  if (coeff.rows() != action_->n_structures() || coeff.cols() <= 0 ||
      !coeff.allFinite()) throw std::invalid_argument("invalid structure block");
  const int na = static_cast<int>(alpha_.size());
  const int nb = static_cast<int>(beta_.size());
  Eigen::MatrixXd result(action_->n_structures(), coeff.cols());
  if (complete_carrier_) {
    ExteriorOrbitalTransform la(alpha_, left_upper), lb(beta_, left_upper);
    ExteriorOrbitalTransform ra(alpha_, right_upper), rb(beta_, right_upper);
    const auto dla = la.direction(delta_left_upper);
    const auto dlb = lb.direction(delta_left_upper);
    const auto dra = ra.direction(delta_right_upper);
    const auto drb = rb.direction(delta_right_upper);
    for (int col = 0; col < coeff.cols(); ++col) {
      Eigen::MatrixXd one = coeff.col(col);
      Eigen::MatrixXd dy;
      Eigen::MatrixXd y = carrier_image(*action_, ra, rb, one,
                                        &dra, &drb, &dy);
      // Differentiate F_left^T y: its seed derivative is dy.
      Eigen::MatrixXd seed_direction = std::move(dy);
      carrier_adjoint(*action_, la, lb,
          std::move(y), &dla, &dlb, &seed_direction);
      result.col(col) = action_->contract_spin_product_block(seed_direction);
    }
    return result;
  }
  const Eigen::MatrixXd k = left_upper.transpose() * right_upper;
  const Eigen::MatrixXd dk = delta_left_upper.transpose() * right_upper +
                             left_upper.transpose() * delta_right_upper;
  const std::size_t denominator = sizeof(double) *
      static_cast<std::size_t>(std::max(1, 4 * nb + 1));
  const int rows_per_tile = std::max(1, std::min(
      na, static_cast<int>(kRowWorkspaceBytes / denominator)));
  for (int col = 0; col < coeff.cols(); ++col) {
    Eigen::MatrixXd one_coeff = coeff.col(col);
    Eigen::MatrixXd one_result = Eigen::MatrixXd::Zero(action_->n_structures(), 1);
    for (int begin = 0; begin < na; begin += rows_per_tile) {
      const int rows = std::min(rows_per_tile, na - begin);
      Eigen::MatrixXd projected = Eigen::MatrixXd::Zero(rows, nb);
      Eigen::MatrixXd dprojected = Eigen::MatrixXd::Zero(rows, nb);
      for (int ar = 0; ar < na; ++ar) {
        const Eigen::MatrixXd source = row_of(*action_, one_coeff, ar, nb);
        for (int ai = 0; ai < rows; ++ai) {
          const auto& left = alpha_[begin + ai];
          const auto& right = alpha_[ar];
          projected.row(ai).noalias() += pair_value(k, left, right) * source.row(0);
          dprojected.row(ai).noalias() +=
              pair_direction(k, dk, left, right) * source.row(0);
        }
      }
      Eigen::MatrixXd tile(rows, nb);
      for (int bl = 0; bl < nb; ++bl) {
        Eigen::VectorXd overlap(nb), doverlap(nb);
        for (int br = 0; br < nb; ++br) {
          overlap[br] = pair_value(k, beta_[bl], beta_[br]);
          doverlap[br] = pair_direction(k, dk, beta_[bl], beta_[br]);
        }
        tile.col(bl).noalias() = dprojected * overlap + projected * doverlap;
      }
      action_->add_spin_product_tile(tile, begin, 0, &one_result);
    }
    result.col(col) = one_result;
  }
  return result;
}

UpperPair RepresentationOverlap::adjoint(
    const Eigen::Ref<const Eigen::MatrixXd>& left_upper,
    const Eigen::Ref<const Eigen::MatrixXd>& right_upper,
    const Eigen::Ref<const Eigen::VectorXd>& left_coeff,
    const Eigen::Ref<const Eigen::VectorXd>& right_coeff) const {
  validate_chart(left_upper, n_active_);
  validate_chart(right_upper, n_active_);
  if (left_coeff.size() != action_->n_structures() ||
      right_coeff.size() != action_->n_structures() ||
      !left_coeff.allFinite() || !right_coeff.allFinite())
    throw std::invalid_argument("invalid structure coefficient vectors");
  const Eigen::MatrixXd left = left_coeff;
  const Eigen::MatrixXd right = right_coeff;
  const int na = static_cast<int>(alpha_.size());
  const int nb = static_cast<int>(beta_.size());
  if (complete_carrier_) {
    ExteriorOrbitalTransform la(alpha_, left_upper), lb(beta_, left_upper);
    ExteriorOrbitalTransform ra(alpha_, right_upper), rb(beta_, right_upper);
    const Eigen::MatrixXd yl = carrier_image(*action_, la, lb, left);
    const Eigen::MatrixXd yr = carrier_image(*action_, ra, rb, right);
    const Eigen::MatrixXd ml = action_->carrier_generator_adjoint(yr, yl);
    const Eigen::MatrixXd mr = action_->carrier_generator_adjoint(yl, yr);
    return {upper_part(ml * chart_inverse_transpose(left_upper)),
            upper_part(mr * chart_inverse_transpose(right_upper))};
  }
  const Eigen::MatrixXd k = left_upper.transpose() * right_upper;
  Eigen::MatrixXd g = Eigen::MatrixXd::Zero(n_active_, n_active_);
  // Project one right-alpha row through the beta overlap, then contract all
  // left-alpha rows. This has no U_alpha^2 U_beta^2 traversal or pair cache.
  for (int ar = 0; ar < na; ++ar) {
    const Eigen::MatrixXd xr = row_of(*action_, right, ar, nb);
    Eigen::VectorXd projected = Eigen::VectorXd::Zero(nb);
    for (int bl = 0; bl < nb; ++bl)
      for (int br = 0; br < nb; ++br)
        projected[bl] += pair_value(k, beta_[bl], beta_[br]) * xr(0, br);
    for (int al = 0; al < na; ++al) {
      const Eigen::MatrixXd xl = row_of(*action_, left, al, nb);
      const double weight = xl.row(0).dot(projected);
      add_pair_gradient(k, nullptr, alpha_[al], alpha_[ar],
                        weight, 0.0, &g, nullptr);
    }
  }
  // The symmetric projection by alpha yields all beta-pair weights.
  for (int br = 0; br < nb; ++br) {
    const Eigen::MatrixXd xr = column_of(*action_, right, br, na);
    Eigen::VectorXd projected = Eigen::VectorXd::Zero(na);
    for (int al = 0; al < na; ++al)
      for (int ar = 0; ar < na; ++ar)
        projected[al] += pair_value(k, alpha_[al], alpha_[ar]) * xr(ar, 0);
    for (int bl = 0; bl < nb; ++bl) {
      const Eigen::MatrixXd xl = column_of(*action_, left, bl, na);
      const double weight = xl.col(0).dot(projected);
      add_pair_gradient(k, nullptr, beta_[bl], beta_[br],
                        weight, 0.0, &g, nullptr);
    }
  }
  return {upper_part(right_upper * g.transpose()),
          upper_part(left_upper * g)};
}

UpperPair RepresentationOverlap::adjoint_direction(
    const Eigen::Ref<const Eigen::MatrixXd>& left_upper,
    const Eigen::Ref<const Eigen::MatrixXd>& right_upper,
    const Eigen::Ref<const Eigen::MatrixXd>& delta_left_upper,
    const Eigen::Ref<const Eigen::MatrixXd>& delta_right_upper,
    const Eigen::Ref<const Eigen::VectorXd>& left_coeff,
    const Eigen::Ref<const Eigen::VectorXd>& right_coeff) const {
  validate_chart(left_upper, n_active_);
  validate_chart(right_upper, n_active_);
  validate_direction(delta_left_upper, n_active_);
  validate_direction(delta_right_upper, n_active_);
  if (left_coeff.size() != action_->n_structures() ||
      right_coeff.size() != action_->n_structures() ||
      !left_coeff.allFinite() || !right_coeff.allFinite())
    throw std::invalid_argument("invalid structure coefficient vectors");
  const Eigen::MatrixXd left = left_coeff;
  const Eigen::MatrixXd right = right_coeff;
  const int na = static_cast<int>(alpha_.size());
  const int nb = static_cast<int>(beta_.size());
  if (complete_carrier_) {
    ExteriorOrbitalTransform la(alpha_, left_upper), lb(beta_, left_upper);
    ExteriorOrbitalTransform ra(alpha_, right_upper), rb(beta_, right_upper);
    const auto dla = la.direction(delta_left_upper);
    const auto dlb = lb.direction(delta_left_upper);
    const auto dra = ra.direction(delta_right_upper);
    const auto drb = rb.direction(delta_right_upper);
    Eigen::MatrixXd dyl, dyr;
    const Eigen::MatrixXd yl = carrier_image(*action_, la, lb, left,
                                             &dla, &dlb, &dyl);
    const Eigen::MatrixXd yr = carrier_image(*action_, ra, rb, right,
                                             &dra, &drb, &dyr);
    const Eigen::MatrixXd ml = action_->carrier_generator_adjoint(yr, yl);
    const Eigen::MatrixXd mr = action_->carrier_generator_adjoint(yl, yr);
    const Eigen::MatrixXd dml =
        action_->carrier_generator_adjoint(dyr, yl) +
        action_->carrier_generator_adjoint(yr, dyl);
    const Eigen::MatrixXd dmr =
        action_->carrier_generator_adjoint(dyl, yr) +
        action_->carrier_generator_adjoint(yl, dyr);
    const Eigen::MatrixXd il = chart_inverse_transpose(left_upper);
    const Eigen::MatrixXd ir = chart_inverse_transpose(right_upper);
    return {upper_part(dml * il - ml * il * delta_left_upper.transpose() * il),
            upper_part(dmr * ir - mr * ir * delta_right_upper.transpose() * ir)};
  }
  const Eigen::MatrixXd k = left_upper.transpose() * right_upper;
  const Eigen::MatrixXd dk = delta_left_upper.transpose() * right_upper +
                             left_upper.transpose() * delta_right_upper;
  Eigen::MatrixXd g = Eigen::MatrixXd::Zero(n_active_, n_active_);
  Eigen::MatrixXd dg = Eigen::MatrixXd::Zero(n_active_, n_active_);
  for (int ar = 0; ar < na; ++ar) {
    const Eigen::MatrixXd xr = row_of(*action_, right, ar, nb);
    Eigen::VectorXd projected = Eigen::VectorXd::Zero(nb);
    Eigen::VectorXd dprojected = Eigen::VectorXd::Zero(nb);
    for (int bl = 0; bl < nb; ++bl)
      for (int br = 0; br < nb; ++br) {
        projected[bl] += pair_value(k, beta_[bl], beta_[br]) * xr(0, br);
        dprojected[bl] += pair_direction(k, dk, beta_[bl], beta_[br]) * xr(0, br);
      }
    for (int al = 0; al < na; ++al) {
      const Eigen::MatrixXd xl = row_of(*action_, left, al, nb);
      add_pair_gradient(k, &dk, alpha_[al], alpha_[ar],
                        xl.row(0).dot(projected),
                        xl.row(0).dot(dprojected), &g, &dg);
    }
  }
  for (int br = 0; br < nb; ++br) {
    const Eigen::MatrixXd xr = column_of(*action_, right, br, na);
    Eigen::VectorXd projected = Eigen::VectorXd::Zero(na);
    Eigen::VectorXd dprojected = Eigen::VectorXd::Zero(na);
    for (int al = 0; al < na; ++al)
      for (int ar = 0; ar < na; ++ar) {
        projected[al] += pair_value(k, alpha_[al], alpha_[ar]) * xr(ar, 0);
        dprojected[al] += pair_direction(k, dk, alpha_[al], alpha_[ar]) * xr(ar, 0);
      }
    for (int bl = 0; bl < nb; ++bl) {
      const Eigen::MatrixXd xl = column_of(*action_, left, bl, na);
      add_pair_gradient(k, &dk, beta_[bl], beta_[br],
                        xl.col(0).dot(projected),
                        xl.col(0).dot(dprojected), &g, &dg);
    }
  }
  return {upper_part(delta_right_upper * g.transpose() +
                     right_upper * dg.transpose()),
          upper_part(delta_left_upper * g + left_upper * dg)};
}

}  // namespace xmvb::vb
