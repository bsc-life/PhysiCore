#pragma once

#include "agent_data.h"
#include "base_agent.h"
#include "types.h"

namespace physicore::biofvm {

template <typename AgentDataType>
class agent_generic : public physicore::base_agent
{
protected:
	AgentDataType& data;

public:
	agent_generic(index_t index, AgentDataType& data) : physicore::base_agent(index, data.base_data), data(data) {}

	std::span<real_t> secretion_rates()
	{
		return std::span<real_t>(&data.secretion_rates[index * data.substrate_count], data.substrate_count);
	}

	std::span<real_t> saturation_densities()
	{
		return std::span<real_t>(&data.saturation_densities[index * data.substrate_count], data.substrate_count);
	}

	std::span<real_t> uptake_rates()
	{
		return std::span<real_t>(&data.uptake_rates[index * data.substrate_count], data.substrate_count);
	}

	std::span<real_t> net_export_rates()
	{
		return std::span<real_t>(&data.net_export_rates[index * data.substrate_count], data.substrate_count);
	}

	std::span<real_t> internalized_substrates()
	{
		return std::span<real_t>(&data.internalized_substrates[index * data.substrate_count], data.substrate_count);
	}

	std::span<real_t> fraction_released_at_death()
	{
		return std::span<real_t>(&data.fraction_released_at_death[index * data.substrate_count], data.substrate_count);
	}

	std::span<real_t> fraction_transferred_when_ingested()
	{
		return std::span<real_t>(&data.fraction_transferred_when_ingested[index * data.substrate_count],
								 data.substrate_count);
	}

	real_t& volume() { return data.volumes[index]; }
};

using agent = agent_generic<agent_data>;

} // namespace physicore::biofvm
