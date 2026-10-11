#pragma once

namespace xmvb::vb {
class StructureAction;
struct AcceptedPointContext;
}

namespace xmvb::tools {

/** @brief Bounded reference audit of response spectra and factorized metric inverses.
 * @param limit Explicit maximum structure and spin-product dimensions.
 * @param complete_space Whether the chosen structures span the complete spin sector.
 * Dense reference allocations are diagnostic only and rejected beyond the limit.
 */
void audit_structure_response(
    const vb::StructureAction& action,
    const vb::AcceptedPointContext& context,
    bool complete_space,
    int limit);

}  // namespace xmvb::tools
