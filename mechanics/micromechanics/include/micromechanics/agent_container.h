#pragma once

#include <memory>

#include <common/base_agent.h>
#include <common/generic_agent_container.h>

#include "agent.h"

namespace physicore::mechanics::micromechanics {

class agent_container : public generic_agent_and_data_container<base_agent, agent>
{
public:
	using generic_agent_and_data_container<base_agent, agent>::generic_agent_and_data_container;
};

} // namespace physicore::mechanics::micromechanics
