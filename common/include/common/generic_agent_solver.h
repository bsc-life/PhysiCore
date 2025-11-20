#pragma once

#include "concepts.h"

namespace physicore {

template <derived_from_base_agent T>
class generic_agent_interface_container;

template <derived_from_base_agent T>
class generic_agent_impl_container;

template <derived_from_base_agent AgentType>
class generic_agent_solver
{
public:
	typename AgentType::DataType& retrieve_agent_data(
		generic_agent_interface_container<typename AgentType::InterfaceType>& container)
	{
		auto& casted_container = dynamic_cast<generic_agent_impl_container<AgentType>&>(container);
		return casted_container.data;
	}
};

} // namespace physicore
