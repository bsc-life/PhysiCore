#pragma once

#include <stdexcept>

#include "agent_data.h"
#include "base_agent.h"
#include "types.h"

namespace physicore::biofvm {

class agent : public physicore::base_agent
{
protected:
	agent_data& data;

public:
	agent(index_t index, agent_data& data);

	std::span<real_t> internalized_substrates();
	std::span<real_t> fraction_released_at_death();
	std::span<real_t> fraction_transferred_when_ingested();

	real_t& volume();
};

template <typename AgentType>
class secretion_uptake_agent_mixin : public AgentType
{
public:
	template <typename... Args>
	secretion_uptake_agent_mixin(Args&&... args) : AgentType(std::forward<Args>(args)...)
	{
		if (this->data.template get<secretion_uptake_component>() == nullptr)
			throw std::runtime_error("Agent data does not contain the required component: secretion_uptake_component");
	}

	std::span<real_t> secretion_rates()
	{
		auto& data = this->data.template get<secretion_uptake_component>();
		return std::span<real_t>(&data.secretion_rates[index * data.substrate_count], data.substrate_count);
	}

	std::span<real_t> saturation_densities()
	{
		auto& data = this->data.template get<secretion_uptake_component>();
		return std::span<real_t>(&data.saturation_densities[index * data.substrate_count], data.substrate_count);
	}

	std::span<real_t> uptake_rates()
	{
		auto& data = this->data.template get<secretion_uptake_component>();
		return std::span<real_t>(&data.uptake_rates[index * data.substrate_count], data.substrate_count);
	}
};

template <typename AgentType>
class net_export_agent_mixin : public AgentType
{
public:
	template <typename... Args>
	net_export_agent_mixin(Args&&... args) : AgentType(std::forward<Args>(args)...)
	{
		if (this->data.template get<net_export_component>() == nullptr)
			throw std::runtime_error("Agent data does not contain the required component: net_export_component");
	}

	std::span<real_t> net_export_rates()
	{
		auto& data = this->data.template get<net_export_component>();
		return std::span<real_t>(&data.net_export_rates[index * data.substrate_count], data.substrate_count);
	}
};

} // namespace physicore::biofvm
