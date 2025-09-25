#include "agent.h"

using namespace physicore;
using namespace physicore::biofvm;

agent::agent(index_t index, agent_data& data) : physicore::base_agent(index, data.base_data), data(data) {}

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
