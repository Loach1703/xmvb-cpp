#pragma once

#include "libcint/materialized_input.hpp"
#include "vbscf/integrals/ao/libcint/input.hpp"

namespace xmvb::vb {

struct LibcintMaterializedIntegralProviderOptions {
  double integral_tolerance = 1.0e-10;
  bool use_historical_shell_pair_prescreen = true;
  /** @brief Disjoint canonical shell-quartet ownership, selected before libcint. */
  int shard_index = 0;
  int shard_count = 1;
  /** @brief Materialize the core Hamiltonian on shard zero only, when requested. */
  bool include_core_hamiltonian = true;
};

/**
 * @brief Generates complete or disjoint canonical AO ERIs with libcint.
 *
 * Quartet ownership precedes evaluation and allocation. A sharded result
 * must be consumed by distributed contractions, not a complete-input loader.
 */
class LibcintMaterializedIntegralProvider {
public:
  MaterializedAoIntegralBuffers build(
      const LibcintInput& input,
      const LibcintMaterializedIntegralProviderOptions& options = {}) const;
};

}  // namespace xmvb::vb
