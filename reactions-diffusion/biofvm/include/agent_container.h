#pragma once

#include "agent.h"
#include "base_agent.h"
#include "generic_agent_container.h"

namespace physicore::biofvm {

using agent_container = physicore::generic_agent_and_data_container<base_agent, agent>;
using agent_container_interface = physicore::generic_agent_interface_container<agent_interface>;

} // namespace physicore::biofvm
