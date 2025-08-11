#pragma once

#include <memory>

#include "base_agent_data.h"

namespace physicore {

class base_agent;

class base_agent_container
{
	base_agent_data data;

	std::vector<std::unique_ptr<base_agent>> agents;

public:
	virtual base_agent* create();

	virtual void remove(base_agent* agent);

	virtual void remove_at(index_t index);

	virtual ~base_agent_container() = default;
};

} // namespace physicore
