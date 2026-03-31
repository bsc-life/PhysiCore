#pragma once

#include <common/base_agent.h>
#include <common/generic_agent_container.h>

#include "agent.h"

namespace physicore::mechanics::micromechanics {

using agent_container = generic_agent_and_data_container<base_agent, agent>;

} // namespace physicore::mechanics::micromechanics
