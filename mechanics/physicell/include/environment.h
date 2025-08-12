#pragma once

#include "process.h"
#include "types.h"

namespace physicore::micromechanics {

class environment : public process
{
	real_t timestep;

public:
	explicit environment(real_t timestep) : timestep(timestep) {}

	void run_single_timestep() override;
};

} // namespace physicore::micromechanics
