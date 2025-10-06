#pragma once

#include "agent_data.h"
#include "agent_generic_storage.h"

namespace physicore::biofvm {

using agent = agent_generic_storage<base_agent_data, agent_data>;

} // namespace physicore::biofvm
