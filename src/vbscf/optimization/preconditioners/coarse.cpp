#include "vbscf/optimization/preconditioners/coarse.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>

namespace xmvb::vb {

CoarsePreconditioner::CoarsePreconditioner(
    int dimension, int capacity, double lindep, double action_accuracy)
    : lindep_(lindep), action_accuracy_(action_accuracy) {
  if (dimension <= 0 || capacity < 0 || !(lindep > 0) ||
      !std::isfinite(action_accuracy) || action_accuracy < 0)
    throw std::invalid_argument("invalid coarse preconditioner dimensions");
  basis_.resize(dimension,std::min(dimension,capacity));
  images_.resize(dimension,basis_.cols());
  metrics_.resize(dimension,basis_.cols());
}

void CoarsePreconditioner::observe(const Eigen::VectorXd& direction,
    const Eigen::VectorXd& image, const Eigen::VectorXd& metric_image) {
  if (direction.size() != basis_.rows() || image.size() != basis_.rows() ||
      metric_image.size() != basis_.rows() || !direction.allFinite() ||
      !image.allFinite() || !metric_image.allFinite())
    throw std::invalid_argument("invalid coarse operator sample");
  if (size_ == basis_.cols()) return;
  Eigen::VectorXd v = direction, h = image, g = metric_image;
  const double initial = v.stableNorm();
  if (!(initial > 0)) return;
  for (int pass = 0; pass < 2; ++pass)
    for (int i = 0; i < size_; ++i) {
      const double coefficient = basis_.col(i).dot(v);
      v.noalias() -= coefficient*basis_.col(i);
      h.noalias() -= coefficient*images_.col(i);
      g.noalias() -= coefficient*metrics_.col(i);
    }
  const double norm = v.stableNorm();
  if (!(norm > std::max(lindep_,64*std::numeric_limits<double>::epsilon())*initial)) return;
  basis_.col(size_) = v/norm;
  images_.col(size_) = h/norm;
  metrics_.col(size_) = g/norm;
  ++size_;
}

Eigen::VectorXd CoarsePreconditioner::apply(const Eigen::VectorXd& residual,
    double eigenvalue,
    const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& inverse) const {
  if (residual.size() != basis_.rows() || !residual.allFinite() || !std::isfinite(eigenvalue))
    throw std::invalid_argument("invalid coarse inverse residual");
  if (!size_) return inverse(residual);
  const auto v = basis_.leftCols(size_);
  const Eigen::MatrixXd lv = images_.leftCols(size_)-eigenvalue*metrics_.leftCols(size_);
  const Eigen::MatrixXd raw = v.transpose()*lv;
  const double skew_bound = 256*std::numeric_limits<double>::epsilon()*size_*
      std::max(1.0,raw.norm())+action_accuracy_*v.norm()*images_.leftCols(size_).norm();
  if ((raw-raw.transpose()).norm() > skew_bound)
    throw std::runtime_error("coarse sampled shifted operator is not symmetric");
  const Eigen::MatrixXd e = 0.5*(raw+raw.transpose());
  Eigen::FullPivLU<Eigen::MatrixXd> factor(e);
  if (!factor.isInvertible())
    throw std::runtime_error("singular exact coarse shifted operator");
  const auto solve = [&](const Eigen::VectorXd& rhs) -> Eigen::VectorXd {
    const Eigen::VectorXd solution = factor.solve(rhs);
    const double bound = 256*std::numeric_limits<double>::epsilon()*size_*
        std::max(1.0,e.norm()*solution.norm()+rhs.norm());
    if (!solution.allFinite() || (e*solution-rhs).norm() > bound)
      throw std::runtime_error("coarse inverse failed its residual certificate");
    return solution;
  };
  const Eigen::VectorXd a = solve(v.transpose()*residual);
  const Eigen::VectorXd fine = inverse(residual-lv*a);
  const Eigen::VectorXd b = solve(lv.transpose()*fine);
  return fine+v*(a-b);
}

CoarseSpectrum CoarsePreconditioner::orbital_spectrum(int no,
    const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& inverse) const {
  if (no <= 0 || no > basis_.rows())
    throw std::invalid_argument("invalid sampled orbital dimension");
  std::vector<int> indices;
  for (int i=0;i<size_;++i)
    if (basis_.col(i).head(no).squaredNorm() > 0.0 &&
        basis_.col(i).tail(basis_.rows()-no).squaredNorm() == 0.0)
      indices.push_back(i);
  CoarseSpectrum result;
  result.dimension = static_cast<int>(indices.size());
  if (!result.dimension) return result;
  Eigen::MatrixXd v(no,result.dimension), av(no,result.dimension), pav(no,result.dimension);
  for (int j=0;j<result.dimension;++j) {
    v.col(j) = basis_.col(indices[j]).head(no);
    av.col(j) = images_.col(indices[j]).head(no);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(basis_.rows());
    rhs.head(no) = av.col(j);
    const Eigen::VectorXd image = inverse(rhs);
    if (image.size() != basis_.rows() || !image.allFinite())
      throw std::runtime_error("invalid sampled base inverse image");
    pav.col(j) = image.head(no);
  }
  Eigen::MatrixXd e = v.transpose()*av;
  e = (0.5*(e+e.transpose())).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> curvature(e,Eigen::EigenvaluesOnly);
  if (curvature.info() != Eigen::Success)
    throw std::runtime_error("sampled orbital spectrum failed");
  result.curvature_min = curvature.eigenvalues()[0];
  result.curvature_max = curvature.eigenvalues().tail(1)[0];
  result.positive_curvature = result.curvature_min > lindep_*result.curvature_max;
  if (!result.positive_curvature) return result;
  Eigen::MatrixXd t = av.transpose()*pav;
  t = (0.5*(t+t.transpose())).eval();
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(t,e,Eigen::EigenvaluesOnly);
  if (spectrum.info() != Eigen::Success)
    throw std::runtime_error("sampled preconditioned spectrum failed");
  result.preconditioned_min = spectrum.eigenvalues()[0];
  result.preconditioned_max = spectrum.eigenvalues().tail(1)[0];
  return result;
}

}  // namespace xmvb::vb
