#pragma once

#include <cassert>
#include <memory>
#include <vector>

#include "base_agent_interface.h"
#include "concepts.h"
#include "types.h"

namespace physicore {

template <derived_from_base_agent AgentType>
class generic_agent_interface_container
{
protected:
	index_t& get_agent_index(base_agent_interface* agent) { return agent->index; }

public:
	virtual AgentType* create() = 0;

	virtual AgentType* get_agent_at(index_t position) = 0;

	virtual void remove_agent(base_agent_interface* agent) = 0;

	virtual void remove_at(index_t position) = 0;

	virtual std::size_t size() const = 0;

	virtual ~generic_agent_interface_container() = default;
};

template <derived_from_base_agent AgentType>
class generic_agent_solver;

template <derived_from_base_agent AgentType>
class generic_agent_impl_container : public generic_agent_interface_container<typename AgentType::InterfaceType>
{
protected:
	friend generic_agent_solver<AgentType>;

	typename AgentType::DataType& data;

public:
	generic_agent_impl_container(AgentType::DataType& data) : data(data) {}
};

template <derived_from_base_agent... AgentTypes>
class generic_agent_and_data_container : public generic_agent_impl_container<AgentTypes>...
{
protected:
	using MostConcreteAgentType = std::tuple_element_t<sizeof...(AgentTypes) - 1, std::tuple<AgentTypes...>>;

	std::vector<std::unique_ptr<MostConcreteAgentType>> agents;

public:
	std::tuple<std::unique_ptr<typename AgentTypes::DataType>...> agent_datas;

	explicit generic_agent_and_data_container(std::unique_ptr<typename AgentTypes::DataType>&&... datas)
		: generic_agent_impl_container<AgentTypes>(*datas)..., agent_datas(std::move(datas)...)
	{}

	MostConcreteAgentType* create() override
	{
		(std::get<std::unique_ptr<typename AgentTypes::DataType>>(agent_datas)->add(), ...);
		auto new_agent = std::make_unique<MostConcreteAgentType>(this->size(), agent_datas);
		MostConcreteAgentType* agent_ptr = new_agent.get();
		agents.emplace_back(std::move(new_agent));
		return agent_ptr;
	}

	void remove_agent(base_agent_interface* agent) override
	{
		index_t index = generic_agent_interface_container<base_agent_interface>::get_agent_index(agent);
		remove_at(index);
	}

	void remove_at(index_t position) override
	{
		assert(position < this->size());

		if (static_cast<size_t>(position) >= agents.size())
			return;

		(std::get<std::unique_ptr<typename AgentTypes::DataType>>(agent_datas)->remove_at(position), ...);

		index_t& index = generic_agent_interface_container<base_agent_interface>::get_agent_index(agents.back().get());
		index = position;
		std::swap(agents[position], agents.back());
		agents.resize(this->size() - 1);
	}

	MostConcreteAgentType* get_agent_at(index_t position) override
	{
		assert(position < agents.size());
		if (position >= agents.size())
			return nullptr;

		return agents[position].get();
	}

	std::size_t size() const override { return agents.size(); }
};

} // namespace physicore
