#include "base_agent_container.h"

#include "base_agent.h"

namespace physicore {

base_agent* base_agent_container::create()
{
	data.add();
	agents.push_back(std::make_unique<base_agent>(data.agents_count, data));
	return agents.back().get();
}

void base_agent_container::remove(base_agent* agent) { remove_at(agent->index); }

void base_agent_container::remove_at(index_t index)
{
	if (static_cast<size_t>(index) >= agents.size())
		return;

	data.remove_at(index);

	agents.back()->index = index;
	std::swap(agents[index], agents.back());
	agents.resize(data.agents_count);
}

} // namespace physicore
