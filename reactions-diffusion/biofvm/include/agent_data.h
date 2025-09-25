#pragma once

#include <vector>

#include "base_agent_data.h"
#include "componented_data.h"
#include "data_component.h"
#include "types.h"

namespace physicore::biofvm {

struct agent_data : componented_data
{
	physicore::base_agent_data& base_data;

	// n * substrate_count
	std::vector<real_t> internalized_substrates;
	std::vector<real_t> fraction_released_at_death;
	std::vector<real_t> fraction_transferred_when_ingested;

	// n
	std::vector<real_t> volumes;

	index_t substrate_count;

	explicit agent_data(physicore::base_agent_data& base_data, index_t substrate_count = 1);

	void add();
	void remove_at(index_t position);
};

struct secretion_uptake_component : data_component
{
	// n * substrate_count
	std::vector<real_t> secretion_rates;
	std::vector<real_t> saturation_densities;
	std::vector<real_t> uptake_rates;

	index_t substrate_count;

	explicit secretion_uptake_component(index_t substrate_count = 1);

	void resize(index_t n) override;
	void remove_at(index_t position) override;
};

struct net_export_component : data_component
{
	// n * substrate_count
	std::vector<real_t> net_export_rates;

	index_t substrate_count;

	explicit net_export_component(index_t substrate_count = 1);

	void resize(index_t n) override;
	void remove_at(index_t position) override;
};

} // namespace physicore::biofvm
