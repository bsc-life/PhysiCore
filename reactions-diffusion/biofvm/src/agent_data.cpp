#include "agent_data.h"

#include <cassert>

using namespace physicore::biofvm;

agent_data::agent_data(physicore::base_agent_data& base_data, index_t substrate_count)
	: base_data(base_data), substrate_count(substrate_count)
{}

void agent_data::add()
{
	base_data.add();

	secretion_rates.resize(base_data.agents_count * substrate_count);
	saturation_densities.resize(base_data.agents_count * substrate_count);
	uptake_rates.resize(base_data.agents_count * substrate_count);
	net_export_rates.resize(base_data.agents_count * substrate_count);

	internalized_substrates.resize(base_data.agents_count * substrate_count);
	fraction_released_at_death.resize(base_data.agents_count * substrate_count);
	fraction_transferred_when_ingested.resize(base_data.agents_count * substrate_count);

	volumes.resize(base_data.agents_count);
}

void agent_data::remove_at(index_t position)
{
	assert(position < base_data.agents_count);

	if (position >= base_data.agents_count)
		return;

	base_data.remove_at(position);

	if (position < base_data.agents_count)
	{
		base_agent_data::move_vector(&secretion_rates[position * substrate_count],
									 &secretion_rates[base_data.agents_count * substrate_count], substrate_count);
		base_agent_data::move_vector(&saturation_densities[position * substrate_count],
									 &saturation_densities[base_data.agents_count * substrate_count], substrate_count);
		base_agent_data::move_vector(&uptake_rates[position * substrate_count],
									 &uptake_rates[base_data.agents_count * substrate_count], substrate_count);
		base_agent_data::move_vector(&net_export_rates[position * substrate_count],
									 &net_export_rates[base_data.agents_count * substrate_count], substrate_count);

		base_agent_data::move_vector(&internalized_substrates[position * substrate_count],
									 &internalized_substrates[base_data.agents_count * substrate_count],
									 substrate_count);
		base_agent_data::move_vector(&fraction_released_at_death[position * substrate_count],
									 &fraction_released_at_death[base_data.agents_count * substrate_count],
									 substrate_count);
		base_agent_data::move_vector(&fraction_transferred_when_ingested[position * substrate_count],
									 &fraction_transferred_when_ingested[base_data.agents_count * substrate_count],
									 substrate_count);

		base_agent_data::move_scalar(&volumes[position], &volumes[base_data.agents_count]);
	}

	secretion_rates.resize(base_data.agents_count * substrate_count);
	saturation_densities.resize(base_data.agents_count * substrate_count);
	uptake_rates.resize(base_data.agents_count * substrate_count);
	net_export_rates.resize(base_data.agents_count * substrate_count);

	internalized_substrates.resize(base_data.agents_count * substrate_count);
	fraction_released_at_death.resize(base_data.agents_count * substrate_count);
	fraction_transferred_when_ingested.resize(base_data.agents_count * substrate_count);

	volumes.resize(base_data.agents_count);
}
