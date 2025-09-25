#pragma once

#include "generic_agent_container.h"

namespace physicore {

template <agent_data_type AgentDataType>
class generic_agent_solver
{
public:
	template <derived_from_base_agent AgentType>
	AgentDataType& retrieve_agent_data(generic_agent_container<AgentDataType, AgentType>& container)
	{
		return container.data;
	}
};

} // namespace physicore
