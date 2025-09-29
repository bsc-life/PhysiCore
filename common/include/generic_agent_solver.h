#pragma once

#include <memory>

#include "generic_agent_container.h"

namespace physicore {

template <agent_data_type AgentDataType>
class generic_agent_solver
{
public:
	template <derived_from_base_agent AgentType>
	AgentDataType& retrieve_agent_data(generic_agent_container<AgentType>& container)
	{
		auto& casted_container = dynamic_cast<generic_agent_and_data_container<AgentDataType, AgentType>&>(container);

		return casted_container.data;
	}
};

} // namespace physicore
