#pragma once

#include <memory>
#include <span>

#include "base_agent_interface.h"
#include "types.h"

namespace physicore {

template <typename AgentDataType>
class base_agent_generic_storage : public virtual base_agent_interface
{
public:
	using DataType = AgentDataType;
	using InterfaceType = base_agent_interface;

	AgentDataType& base_data;

public:
	base_agent_generic_storage(index_t id, AgentDataType& data) : base_agent_interface(id), base_data(data)
	{
		assert(index < base_data.agents_count);
	}

	base_agent_generic_storage(index_t id, std::tuple<std::unique_ptr<AgentDataType>>& data)
		: base_agent_generic_storage(id, *std::get<0>(data))
	{}

	std::span<real_t> position() override
	{
		return std::span<real_t>(&base_data.positions[index * base_data.dims], base_data.dims);
	}
};

} // namespace physicore
