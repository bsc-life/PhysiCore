#pragma once

#include <cassert>
#include <ranges>

#include "base_agent_container.h"

namespace physicore {

template <typename AgentDataType, typename AgentType>
class generic_agent_container : public base_agent_container
{
protected:
	AgentDataType data;

public:
	template <typename... Args>
	generic_agent_container(Args&&... args) : base_agent_container(), data(this->base_data, std::forward<Args>(args)...)
	{}

	AgentType* create()
	{
		data.add();
		auto new_agent = std::make_unique<AgentType>(this->base_data.agents_count - 1, data);
		AgentType* agent_ptr = new_agent.get();
		this->agents.emplace_back(std::move(new_agent));
		return agent_ptr;
	}

	AgentType* get_agent_at(index_t position)
	{
		base_agent* base_agent = base_agent_container::get_agent_at(position);
		AgentType* agent = dynamic_cast<AgentType*>(base_agent);

		assert(base_agent == nullptr || agent != nullptr); // if base_agent is not null, dynamic_cast must succeed

		return agent;
	}

	virtual void remove_at(index_t position) override
	{
		assert(position < this->base_data.agents_count);

		if (static_cast<size_t>(position) >= this->agents.size())
			return;

		data.remove(position);

		swap_and_erase_agent(position);
	}

	auto get_agents() const
	{
		return this->agents | std::views::transform([](auto& ptr) { return dynamic_cast<AgentType*>(ptr.get()); });
	}
};

} // namespace physicore
