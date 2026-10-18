#include "integrals/active/raw.hpp"

#include <stdexcept>
#include <vector>

#include "integrals/active/ri_transform.hpp"
#include "integrals/active/transform.hpp"

namespace xmvb::integrals {
namespace {

Eigen::MatrixXd unpack_kernel(const std::vector<double>& packed,
                              Eigen::Index pairs) {
  Eigen::MatrixXd kernel(pairs, pairs);
  for (Eigen::Index p = 0; p < pairs; ++p) {
    for (Eigen::Index q = 0; q <= p; ++q) {
      const double value = packed[static_cast<std::size_t>(p * (p + 1) / 2 + q)];
      kernel(p, q) = value;
      kernel(q, p) = value;
    }
  }
  return kernel;
}

std::vector<double> pack_kernel(const Eigen::MatrixXd& kernel) {
  const Eigen::Index pairs = kernel.rows();
  std::vector<double> packed(
      static_cast<std::size_t>(pairs * (pairs + 1) / 2));
  for (Eigen::Index p = 0; p < pairs; ++p) {
    for (Eigen::Index q = 0; q <= p; ++q) {
      packed[static_cast<std::size_t>(p * (p + 1) / 2 + q)] = kernel(p, q);
    }
  }
  return packed;
}

}  // namespace

std::vector<double> materialize_raw_pair_kernel(
    const std::vector<double>& packed,
    const Eigen::MatrixXd& upper) {
  if (upper.size() == 0) {
    return packed;
  }
  const Eigen::Index n = upper.rows();
  if (n <= 0 || upper.cols() != n) {
    throw std::invalid_argument("raw pair-kernel frame must be square");
  }
  const Eigen::Index pairs = n * (n + 1) / 2;
  if (packed.size() != static_cast<std::size_t>(pairs * (pairs + 1) / 2)) {
    throw std::invalid_argument("raw pair-kernel tensor size is inconsistent");
  }
  return pack_kernel(transform_pair_kernel(unpack_kernel(packed, pairs), upper));
}

ActiveSpaceTwoElectronResult materialize_raw_active_integrals(
    const ActiveSpaceTwoElectronResult& framed) {
  if (framed.orbital_transform.size() == 0) {
    return framed;
  }
  const Eigen::MatrixXd& r = framed.orbital_transform;
  const Eigen::Index n = r.rows();
  if (n <= 0 || r.cols() != n ||
      framed.dense_active_coefficients.rows() <= 0 ||
      framed.dense_active_coefficients.cols() != n) {
    throw std::invalid_argument("active integral frame dimensions are inconsistent");
  }
  const Eigen::Index pairs = n * (n + 1) / 2;
  const auto& packed = framed.packed_active_two_electron_integrals;
  if (!packed.empty() && packed.size() !=
          static_cast<std::size_t>(pairs * (pairs + 1) / 2)) {
    throw std::invalid_argument("framed active integral tensor size is inconsistent");
  }
  if (framed.representation == ActiveSpaceTwoElectronRepresentation::PackedExact &&
      packed.empty()) {
    throw std::invalid_argument("framed exact active integral tensor is missing");
  }
  if (framed.representation == ActiveSpaceTwoElectronRepresentation::ResolutionOfIdentity &&
      (framed.n_auxiliary_functions <= 0 ||
       framed.ri_active_pair_factors.rows() != framed.n_auxiliary_functions ||
       framed.ri_active_pair_factors.cols() != pairs)) {
    throw std::invalid_argument("framed RI active pair factors are inconsistent");
  }

  // Construct only persistent data; AO-pair products and one-index RI images
  // belong to the barred frame and can be much larger than the active tensor.
  ActiveSpaceTwoElectronResult raw;
  raw.representation = framed.representation;
  raw.n_auxiliary_functions = framed.n_auxiliary_functions;
  if (!packed.empty()) {
    raw.packed_active_two_electron_integrals =
        materialize_raw_pair_kernel(packed, r);
  }
  if (framed.representation == ActiveSpaceTwoElectronRepresentation::ResolutionOfIdentity) {
    raw.ri_active_pair_factors =
        transform_ri_pair_factors(framed.ri_active_pair_factors, r);
  }
  raw.dense_active_coefficients = framed.dense_active_coefficients * r;
  return raw;
}

}  // namespace xmvb::integrals
