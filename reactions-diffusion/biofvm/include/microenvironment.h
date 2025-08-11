#pragma once

#include <process.h>
#include <types.h>

namespace physicore::biofvm {

class microenvironment : public process
{
	real_t timestep;

public:
	explicit microenvironment(real_t timestep);

	void run_single_timestep() override;
};

} // namespace physicore::biofvm
