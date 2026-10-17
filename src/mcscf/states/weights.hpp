#pragma once

#include <vector>

namespace xmvb::mcscf {

/** @brief Validate a nonempty list of nonnegative state weights. */
void validate_weights(const std::vector<double>& weights);

/** @brief Normalize state weights without changing their ordering.
 * The sum must be positive. Zero-weight states are retained, not discarded.
 */
std::vector<double> normalize_weights(const std::vector<double>& weights);

} // namespace xmvb::mcscf
