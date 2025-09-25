#pragma once

#include <cassert>
#include <concepts>
#include <ranges>

#include "base_agent_container.h"
#include "types.h"

namespace physicore {

// Concept to check if a type has the required agent data methods
template <typename T>
concept agent_data_type = requires(T data, index_t pos) {
	{ data.add() } -> std::same_as<void>;
	{ data.remove_at(pos) } -> std::same_as<void>;
};

// Concept to check if a type is derived from base_agent
template <typename T>
concept derived_from_base_agent = std::derived_from<T, base_agent>;

template <agent_data_type AgentDataType>
class generic_agent_solver;

template <agent_data_type AgentDataType, derived_from_base_agent AgentType>
class generic_agent_container : public base_agent_container
{
protected:
	AgentDataType data;

	friend generic_agent_solver<AgentDataType>;

public:
	template <typename... Args>
	explicit generic_agent_container(index_t dims = 3, Args&&... args)
		: base_agent_container(dims), data(this->base_data, std::forward<Args>(args)...)
	{}

	AgentType* create() override
	{
		data.add();
		auto new_agent = std::make_unique<AgentType>(this->base_data.agents_count - 1, data);
		AgentType* agent_ptr = new_agent.get();
		this->agents.emplace_back(std::move(new_agent));
		return agent_ptr;
	}

	AgentType* get_agent_at(index_t position) override
	{
		base_agent* base_agent = base_agent_container::get_agent_at(position);
		auto agent = dynamic_cast<AgentType*>(base_agent);

		assert(base_agent == nullptr || agent != nullptr); // if base_agent is not null, dynamic_cast must succeed

		return agent;
	}

	void remove_at(index_t position) override
	{
		assert(position < this->base_data.agents_count);

		if (static_cast<size_t>(position) >= this->agents.size())
			return;

		data.remove_at(position);

		swap_and_erase_agent(position);
	}

	auto get_agents() const
	{
		return this->agents | std::views::transform([](auto& ptr) { return dynamic_cast<AgentType*>(ptr.get()); });
	}
};

} // namespace physicore
