#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace xmvb::mcscf {

/** @brief Number of fixed-particle strings; zero for an absent sector. */
std::size_t string_count(int n_orbitals, int n_electrons);

/** @brief Complete fixed-particle spin strings, with orbital zero as bit zero. */
class SpinStrings {
public:
  SpinStrings(int n_orbitals, int n_electrons);
  const std::vector<std::vector<int>>& occupations() const noexcept { return occupations_; }
  const std::vector<std::uint64_t>& masks() const noexcept { return masks_; }
  int index(std::uint64_t mask) const { return indices_.at(mask); }
  int size() const noexcept { return static_cast<int>(masks_.size()); }

private:
  std::vector<std::uint64_t> masks_;
  std::vector<std::vector<int>> occupations_;
  std::unordered_map<std::uint64_t, int> indices_;
};

}  // namespace xmvb::mcscf
