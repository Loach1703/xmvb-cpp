#include "dft/integration/density.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace xmvb::dft {

void evaluate_density_tile(
    const Eigen::Ref<const Eigen::MatrixXd>& values,
    const Eigen::Ref<const Eigen::MatrixXd>& density_matrix,
    Eigen::MatrixXd* projected,
    Eigen::VectorXd* density) {
  if (!projected || !density ||
      density_matrix.rows() != density_matrix.cols() ||
      values.cols() != density_matrix.rows()) {
    throw std::invalid_argument("AO density tile shape mismatch");
  }
  projected->resize(values.rows(), values.cols());
  projected->noalias() = values * density_matrix;
  density->resize(values.rows());
  for (Eigen::Index point = 0; point < values.rows(); ++point) {
    const auto products =
        (values.row(point).array() * projected->row(point).array());
    const double rho = products.sum();
    const double absolute_sum = products.abs().sum();
    const double noise = 64.0 * std::numeric_limits<double>::epsilon() *
                         absolute_sum;
    if (!std::isfinite(rho) || !std::isfinite(absolute_sum)) {
      throw std::runtime_error("nonfinite AO density on DFT grid");
    }
    if (rho < -noise) {
      throw std::invalid_argument("restricted AO density is negative on DFT grid");
    }
    (*density)[point] = std::max(0.0, rho);
  }
}

}  // namespace xmvb::dft
