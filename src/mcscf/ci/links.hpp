#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <Eigen/Core>

namespace xmvb::mcscf {

/** @brief One exact single replacement, in packed-pair order.
 * The current contract has at most63 orbitals, so pair/orbital indices fit
 * uint16/uint8. Keep the numerical sign in double: node-wise kernels consume
 * it directly. String addresses retain int; no floating precision is narrowed.
 */
struct SpinLink {
  int source;
  std::uint16_t pair;
  std::uint8_t created,annihilated;
  double sign;
};

/** @brief Borrowed immutable string row; no per-field maps or allocation. */
struct DensityConnections {
  const SpinLink* data;
  int count;
  int size() const noexcept { return count; }
  const SpinLink& operator[](int i) const noexcept { return data[i]; }
  const SpinLink* begin() const noexcept { return data; }
  const SpinLink* end() const noexcept { return count ? data+count:data; }
};

/** @brief Fixed-stride topology for a complete one-spin sector.
 * Each string has ne*(norb-ne) links. A flat record array stores those rows;
 * the occupied-orbital Eigen array has one contiguous column per string.
 * Sigma directions and identical alpha/beta sectors share this immutable data.
 */
struct SpinConnections {
  int strings=0,degree=0;
  std::vector<SpinLink> links;
  Eigen::Matrix<std::uint8_t,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> occupations;

  int n_strings() const noexcept { return strings; }
  int n_occupied() const noexcept { return occupations.rows(); }
  int n_singles() const noexcept { return degree; }

  DensityConnections singles(int string) const noexcept {
    return {degree ? links.data()+static_cast<std::size_t>(degree)*string:nullptr,degree};
  }
  Eigen::Map<const Eigen::Matrix<std::uint8_t,Eigen::Dynamic,1>> occupied(int string) const noexcept {
    const auto count=occupations.rows();
    return {count ? occupations.data()+count*string:nullptr,count};
  }
  std::size_t dynamic_bytes() const noexcept {
    return links.capacity()*sizeof(SpinLink)+occupations.size()*sizeof(std::uint8_t);
  }
};

} // namespace xmvb::mcscf
