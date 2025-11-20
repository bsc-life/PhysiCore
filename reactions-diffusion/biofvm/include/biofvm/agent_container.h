#pragma once

#include <common/base_agent.h>
#include <common/generic_agent_container.h>

#include "agent.h"

namespace physicore::biofvm {

using agent_container = physicore::generic_agent_and_data_container<base_agent, agent>;
using agent_container_interface = physicore::generic_agent_interface_container<agent_interface>;
using container_ptr = std::shared_ptr<agent_container_interface>;

} // namespace physicore::biofvm
