#pragma once

#include "agent.h"
#include "generic_agent_container.h"

namespace physicore::biofvm {

using agent_container = physicore::generic_agent_and_data_container<agent_data, agent>;
using agent_container_base = physicore::generic_agent_container<agent>;

} // namespace physicore::biofvm
