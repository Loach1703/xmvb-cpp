#pragma once

#include <cstddef>
#include <vector>

namespace xmvb::vb {

/** @brief Exact source-column index of the original touched spin products.
 * Products are sorted by alpha then beta; no orthogonal product is generated.
 */
struct SpinSources {
  std::vector<int> columns, representatives, alpha_rows, alpha_offsets;
  std::vector<std::vector<int>> beta_groups;

  SpinSources(const std::vector<int>& products, int na, int nb, int parity = 0)
      : alpha_offsets(na + 1, 0), beta_groups(nb) {
    std::vector<unsigned char> source(nb, 0);
    for (int group = 0; group < static_cast<int>(products.size()); ++group) {
      const int a = products[group] / nb, b = products[group] % nb;
      ++alpha_offsets[a + 1];
      beta_groups[b].push_back(group);
      if (!parity || a < b || (a == b && parity > 0)) source[b] = 1;
    }
    for (int a = 0; a < na; ++a) {
      if (alpha_offsets[a + 1]) alpha_rows.push_back(a);
      alpha_offsets[a + 1] += alpha_offsets[a];
    }
    for (int b = 0; b < nb; ++b) {
      if (!beta_groups[b].empty()) columns.push_back(b);
      if (source[b]) representatives.push_back(b);
    }
  }

  std::size_t dynamic_bytes() const noexcept {
    std::size_t bytes = (columns.capacity() + representatives.capacity() + alpha_rows.capacity() + alpha_offsets.capacity()) * sizeof(int) +
        beta_groups.capacity() * sizeof(std::vector<int>);
    for (const auto& groups : beta_groups) bytes += groups.capacity() * sizeof(int);
    return bytes;
  }
};

} // namespace xmvb::vb
