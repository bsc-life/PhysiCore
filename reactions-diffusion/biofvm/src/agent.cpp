#include "agent.h"

using namespace physicore;
using namespace physicore::biofvm;

agent::agent(agent_data& data, index_t index) : physicore::base_agent(index, data.base_data), data(data) {}

std::span<real_t> agent::secretion_rates()
{
	return std::span<real_t>(&data.secretion_rates[index * data.substrate_count], data.substrate_count);
}

std::span<real_t> agent::saturation_densities()
{
	return std::span<real_t>(&data.saturation_densities[index * data.substrate_count], data.substrate_count);
}

std::span<real_t> agent::uptake_rates()
{
	return std::span<real_t>(&data.uptake_rates[index * data.substrate_count], data.substrate_count);
}

std::span<real_t> agent::net_export_rates()
{
	return std::span<real_t>(&data.net_export_rates[index * data.substrate_count], data.substrate_count);
}

std::span<real_t> agent::internalized_substrates()
{
	return std::span<real_t>(&data.internalized_substrates[index * data.substrate_count], data.substrate_count);
}

std::span<real_t> agent::fraction_released_at_death()
{
	return std::span<real_t>(&data.fraction_released_at_death[index * data.substrate_count], data.substrate_count);
}

std::span<real_t> agent::fraction_transferred_when_ingested()
{
	return std::span<real_t>(&data.fraction_transferred_when_ingested[index * data.substrate_count],
							 data.substrate_count);
}

real_t& agent::volume() { return data.volumes[index]; }
