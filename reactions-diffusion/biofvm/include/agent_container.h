#pragma once

#include "agent.h"
#include "generic_agent_container.h"

namespace physicore::biofvm {

using agent_container = physicore::generic_agent_container<agent_data, agent>;

} // namespace physicore::biofvm
