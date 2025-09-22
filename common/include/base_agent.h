#pragma once

#include <span>

#include "base_agent_data.h"

namespace physicore {

class base_agent_container;

class base_agent
{
protected:
	index_t index;
	base_agent_data& base_data;

	friend base_agent_container;

public:
	base_agent(index_t index, base_agent_data& data);

	std::span<real_t> position();

	virtual ~base_agent() = default;
};

} // namespace physicore
