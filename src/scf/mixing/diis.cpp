#include "scf/mixing/diis.hpp"

#include <stdexcept>

#include <Eigen/LU>

namespace xmvb::scf {

Eigen::MatrixXd commutator_error(
    const Eigen::Ref<const Eigen::MatrixXd>& fock,
    const Eigen::Ref<const Eigen::MatrixXd>& density,
    const Eigen::Ref<const Eigen::MatrixXd>& overlap) {
  if (fock.rows() != fock.cols() ||
      density.rows() != fock.rows() || density.cols() != fock.cols() ||
      overlap.rows() != fock.rows() || overlap.cols() != fock.cols()) {
    throw std::invalid_argument("SCF commutator matrix dimensions do not match");
  }
  const Eigen::MatrixXd fps = fock * density * overlap;
  return fps - fps.transpose();
}

Diis::Diis(int max_history) : max_history_(max_history) {
  if (max_history < 0) {
    throw std::invalid_argument("SCF DIIS history must be non-negative");
  }
  history_.reserve(max_history);
}

void Diis::push(
    const Eigen::Ref<const Eigen::MatrixXd>& fock,
    const Eigen::Ref<const Eigen::MatrixXd>& error) {
  if (max_history_ == 0) {
    return;
  }
  if (fock.rows() != fock.cols() ||
      error.rows() != fock.rows() || error.cols() != fock.cols()) {
    throw std::invalid_argument("SCF DIIS matrix dimensions do not match");
  }
  if (static_cast<int>(history_.size()) == max_history_) {
    history_.erase(history_.begin());
  }
  history_.push_back({fock, error});
}

bool Diis::extrapolate(Eigen::MatrixXd* fock) const {
  if (fock == nullptr) {
    throw std::invalid_argument("SCF DIIS output must not be null");
  }
  const int size = static_cast<int>(history_.size());
  if (size < 2) {
    return false;
  }

  Eigen::MatrixXd system = Eigen::MatrixXd::Zero(size + 1, size + 1);
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(size + 1);
  rhs(size) = -1.0;
  for (int row = 0; row < size; ++row) {
    for (int column = 0; column < size; ++column) {
      system(row, column) =
          (history_[row].error.array() * history_[column].error.array()).sum();
    }
    system(row, size) = -1.0;
    system(size, row) = -1.0;
  }

  Eigen::FullPivLU<Eigen::MatrixXd> lu(system);
  if (!lu.isInvertible()) {
    return false;
  }
  const Eigen::VectorXd weights = lu.solve(rhs);
  if (!weights.allFinite()) {
    return false;
  }
  *fock = Eigen::MatrixXd::Zero(
      history_.front().fock.rows(), history_.front().fock.cols());
  for (int entry = 0; entry < size; ++entry) {
    *fock += weights(entry) * history_[entry].fock;
  }
  return fock->allFinite();
}

}  // namespace xmvb::scf
