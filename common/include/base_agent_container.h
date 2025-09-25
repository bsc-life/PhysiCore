#pragma once

#include <memory>

#include "base_agent_data.h"

namespace physicore {

class base_agent;

class base_agent_container
{
protected:
	base_agent_data base_data;

	std::vector<std::unique_ptr<base_agent>> agents;

	void swap_and_erase_agent(index_t position);

public:
	explicit base_agent_container(index_t dims = 3);

	virtual base_agent* create();

	virtual base_agent* get_agent_at(index_t position);

	virtual void remove_agent(base_agent* agent);

	virtual void remove_at(index_t position);

	virtual std::size_t size() const;

	virtual ~base_agent_container() = default;
};

} // namespace physicore
