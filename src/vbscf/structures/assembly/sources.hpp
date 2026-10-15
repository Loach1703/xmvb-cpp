#pragma once

#include <cstddef>
#include <vector>

#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"

namespace xmvb::vb {

/** @brief One exact source factor anchored on a row or column of raw support. */
struct SpinSource {
  int string = 0;
  bool beta = true;
};

/** @brief Minimum row/column cover of the original touched spin products.
 * Each input product belongs to exactly one source. The cover minimizes the
 * number of coordinate-anchored outer products, without inspecting amplitudes.
 * Original output support remains separate from an optional parity quotient.
 */
struct SpinSources {
  std::vector<int> columns, alpha_rows, alpha_offsets;
  std::vector<std::vector<int>> beta_groups;
  std::vector<SpinSource> nodes;
  std::vector<int> offsets, groups;
  ExteriorProjection alpha_projection, beta_projection;

  SpinSources(const std::vector<int>& products, int na, int nb, int parity = 0);

  std::size_t dynamic_bytes() const noexcept {
    std::size_t bytes = (columns.capacity() + alpha_rows.capacity() + alpha_offsets.capacity() +
        offsets.capacity() + groups.capacity()) * sizeof(int) +
        beta_groups.capacity() * sizeof(std::vector<int>) + nodes.capacity() * sizeof(SpinSource);
    for (const auto& groups : beta_groups) bytes += groups.capacity() * sizeof(int);
    for (const auto* plan : {&alpha_projection, &beta_projection}) {
      bytes += plan->shear_pairs.capacity() * sizeof(std::vector<int>);
      for (const auto& pairs : plan->shear_pairs) bytes += pairs.capacity() * sizeof(int);
    }
    return bytes;
  }
};

} // namespace xmvb::vb
