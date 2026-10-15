#include "mcscf/ci/strings.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace xmvb::mcscf {

std::size_t string_count(int n, int electrons) {
  if (n<0 || electrons<0 || electrons>n) return 0;
  std::size_t count = 1;
  const int k = std::min(electrons, n - electrons);
  for (int i = 1; i <= k; ++i) {
    if (count > std::numeric_limits<std::size_t>::max() / (n - k + i))
      throw std::length_error("spin-string count overflows");
    count = count * (n - k + i) / i;
  }
  return count;
}

SpinStrings::SpinStrings(int n, int electrons) {
  if (n < 1 || n > 63 || electrons < 0 || electrons > n)
    throw std::invalid_argument("invalid fixed-particle spin space");
  const auto count=string_count(n,electrons);
  if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::length_error("spin-string count exceeds index range");
  masks_.reserve(count);
  occupations_.reserve(count);
  std::uint64_t mask = electrons ? (std::uint64_t{1} << electrons) - 1 : 0;
  const std::uint64_t limit = std::uint64_t{1} << n;
  for (std::size_t i = 0; i < count; ++i) {
    masks_.push_back(mask);
    std::vector<int> occupied;
    for (int orbital = 0; orbital < n; ++orbital)
      if (mask & (std::uint64_t{1} << orbital)) occupied.push_back(orbital);
    occupations_.push_back(std::move(occupied));
    indices_.emplace(mask, static_cast<int>(i));
    if (i + 1 == count) break;
    const std::uint64_t low = mask & -mask;
    const std::uint64_t next = mask + low;
    mask = next | (((next ^ mask) >> 2) / low);
    if (mask >= limit) throw std::logic_error("inconsistent spin-string enumeration");
  }
}

}  // namespace xmvb::mcscf
