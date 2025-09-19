#pragma once

#include "agent_data.h"
#include "base_agent.h"
#include "types.h"

namespace physicore::biofvm {

class agent : public physicore::base_agent
{
protected:
	agent_data& data;

public:
	agent(agent_data& data, index_t index);

	std::span<real_t> secretion_rates();
	std::span<real_t> saturation_densities();
	std::span<real_t> uptake_rates();
	std::span<real_t> net_export_rates();

	std::span<real_t> internalized_substrates();
	std::span<real_t> fraction_released_at_death();
	std::span<real_t> fraction_transferred_when_ingested();

	real_t& volume();
};

} // namespace physicore::biofvm
