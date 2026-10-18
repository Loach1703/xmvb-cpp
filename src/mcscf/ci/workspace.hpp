#pragma once

#include <atomic>
#include <cstddef>
#include <stdexcept>

namespace xmvb::mcscf {

/** @brief Process-wide H/S tile payload budget, configured before numerical work.
 * Input/output residency, sigma and transformation scratch share this budget;
 * it is not an RSS limit. Integrals, immutable links and Davidson vectors are
 * separate. Increasing the budget never changes the selected physical space.
 */
inline constexpr std::size_t kDefaultCiWorkspaceBytes=4ULL*1024*1024*1024;
namespace detail {
inline std::atomic<std::size_t> ci_workspace_limit{kDefaultCiWorkspaceBytes};
}

inline std::size_t ci_workspace_bytes() noexcept {
  return detail::ci_workspace_limit.load(std::memory_order_relaxed);
}

inline void set_ci_workspace_bytes(std::size_t bytes) {
  if (!bytes) throw std::invalid_argument("CI workspace must be positive");
  detail::ci_workspace_limit.store(bytes,std::memory_order_relaxed);
}

} // namespace xmvb::mcscf
