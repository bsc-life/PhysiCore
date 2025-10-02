#pragma once

#include "base_agent_interface.h"
#include "types.h"

namespace physicore::biofvm {

class agent_interface : public virtual base_agent_interface
{
public:
	virtual std::span<real_t> secretion_rates() = 0;

	virtual std::span<real_t> saturation_densities() = 0;

	virtual std::span<real_t> uptake_rates() = 0;

	virtual std::span<real_t> net_export_rates() = 0;

	virtual std::span<real_t> internalized_substrates() = 0;

	virtual std::span<real_t> fraction_released_at_death() = 0;

	virtual std::span<real_t> fraction_transferred_when_ingested() = 0;

	virtual real_t& volume() = 0;
};

} // namespace physicore::biofvm
