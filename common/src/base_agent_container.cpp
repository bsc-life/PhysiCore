#include "base_agent_container.h"

#include <cassert>

#include "base_agent.h"

namespace physicore {

base_agent* base_agent_container::create()
{
	base_data.add();
	agents.push_back(std::make_unique<base_agent>(base_data.agents_count - 1, base_data));
	return agents.back().get();
}

base_agent* base_agent_container::get_agent_at(index_t position)
{
	assert(position < agents.size());
	if (position >= agents.size())
		return nullptr;

	return agents[position].get();
}

void base_agent_container::remove_agent(base_agent* agent) { remove_at(agent->index); }

void base_agent_container::remove_at(index_t position)
{
	assert(position < base_data.agents_count);

	if (static_cast<size_t>(position) >= agents.size())
		return;

	base_data.remove_at(position);

	swap_and_erase_agent(position);
}

void base_agent_container::swap_and_erase_agent(index_t position)
{
	agents.back()->index = position;
	std::swap(agents[position], agents.back());
	agents.resize(base_data.agents_count);
}

std::size_t base_agent_container::size() const { return agents.size(); }

} // namespace physicore
