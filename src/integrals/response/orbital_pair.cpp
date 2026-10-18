#include "integrals/response/orbital_pair.hpp"
#include "core/timing.hpp"

namespace xmvb::integrals {

Eigen::MatrixXd pullback_orbital_pairs(
    const OrbitalPairMatrix& products,
    const OrbitalPairMatrix& pair_adjoint,
    std::span<const int> first,
    std::span<const int> second,
    int n_active) {
  profiling::ScopedPhase phase(profiling::Phase::OrbitalPairPullback);
  const Eigen::Index pairs = first.size();
  Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(n_active * pairs, n_active);
  for (Eigen::Index p = 0; p < pairs; ++p) {
    for (Eigen::Index q = 0; q < pairs; ++q) {
      weights(second[p] * pairs + q, first[p]) += pair_adjoint(q, p);
      weights(first[p] * pairs + q, second[p]) += pair_adjoint(q, p);
    }
  }
  const Eigen::Map<const OrbitalPairMatrix> image(
      products.data(), products.rows() / n_active, n_active * pairs);
  return image * weights;
}

}  // namespace xmvb::integrals
