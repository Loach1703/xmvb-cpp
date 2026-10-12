#pragma once

#include "tnvb/communication/response.hpp"
#include "vbscf/derivatives/gradient/active_space/result.hpp"

namespace xmvb::tnvb {

/** Map P/Pi/D3G to As, Ah and packed Ag, not to a final orbital gradient. */
void populate_active_space_adjoints(
    const WorkerResponse&, vb::ActiveSpaceGradientResult*);

}  // namespace xmvb::tnvb
