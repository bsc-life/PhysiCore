#pragma once

#include <base_agent_data.h>
#include <vector>

#include "types.h"

namespace physicore::biofvm {

template <template <typename...> typename ContainerType = std::vector>
struct agent_data_generic_storage
{
public:
	physicore::base_agent_data_generic_storage<ContainerType>& base_data;

	// n * substrate_count
	ContainerType<real_t> secretion_rates;
	ContainerType<real_t> saturation_densities;
	ContainerType<real_t> uptake_rates;
	ContainerType<real_t> net_export_rates;

	// n * substrate_count
	ContainerType<real_t> internalized_substrates;
	ContainerType<real_t> fraction_released_at_death;
	ContainerType<real_t> fraction_transferred_when_ingested;

	// n
	ContainerType<real_t> volumes;

	index_t agents_count = 0;
	index_t substrate_count;

	explicit agent_data_generic_storage(physicore::base_agent_data_generic_storage<ContainerType>& base_data,
										index_t substrate_count = 1);

	void add();
	void remove_at(index_t position);
};

template <template <typename...> typename ContainerType>
agent_data_generic_storage<ContainerType>::agent_data_generic_storage(
	physicore::base_agent_data_generic_storage<ContainerType>& base_data, index_t substrate_count)
	: base_data(base_data), substrate_count(substrate_count)
{}

template <template <typename...> typename ContainerType>
void agent_data_generic_storage<ContainerType>::add()
{
	++agents_count;

	secretion_rates.resize(agents_count * substrate_count);
	saturation_densities.resize(agents_count * substrate_count);
	uptake_rates.resize(agents_count * substrate_count);
	net_export_rates.resize(agents_count * substrate_count);

	internalized_substrates.resize(agents_count * substrate_count);
	fraction_released_at_death.resize(agents_count * substrate_count);
	fraction_transferred_when_ingested.resize(agents_count * substrate_count);

	volumes.resize(agents_count);
}

template <template <typename...> typename ContainerType>
void agent_data_generic_storage<ContainerType>::remove_at(index_t position)
{
	assert(position < agents_count);

	if (position >= agents_count)
		return;
	--agents_count;

	if (position < agents_count)
	{
		base_agent_data::move_vector(&secretion_rates[position * substrate_count],
									 &secretion_rates[agents_count * substrate_count], substrate_count);
		base_agent_data::move_vector(&saturation_densities[position * substrate_count],
									 &saturation_densities[agents_count * substrate_count], substrate_count);
		base_agent_data::move_vector(&uptake_rates[position * substrate_count],
									 &uptake_rates[agents_count * substrate_count], substrate_count);
		base_agent_data::move_vector(&net_export_rates[position * substrate_count],
									 &net_export_rates[agents_count * substrate_count], substrate_count);

		base_agent_data::move_vector(&internalized_substrates[position * substrate_count],
									 &internalized_substrates[agents_count * substrate_count], substrate_count);
		base_agent_data::move_vector(&fraction_released_at_death[position * substrate_count],
									 &fraction_released_at_death[agents_count * substrate_count], substrate_count);
		base_agent_data::move_vector(&fraction_transferred_when_ingested[position * substrate_count],
									 &fraction_transferred_when_ingested[agents_count * substrate_count],
									 substrate_count);

		base_agent_data::move_scalar(&volumes[position], &volumes[agents_count]);
	}

	secretion_rates.resize(agents_count * substrate_count);
	saturation_densities.resize(agents_count * substrate_count);
	uptake_rates.resize(agents_count * substrate_count);
	net_export_rates.resize(agents_count * substrate_count);

	internalized_substrates.resize(agents_count * substrate_count);
	fraction_released_at_death.resize(agents_count * substrate_count);
	fraction_transferred_when_ingested.resize(agents_count * substrate_count);

	volumes.resize(agents_count);
}

} // namespace physicore::biofvm
