#pragma once

#include <process.h>
#include <types.h>

namespace physicore {
namespace biofvm {

class microenvironment : public process
{
	real_t timestep;

public:
	microenvironment(real_t timestep);

	void run_single_timestep() override;
};

} // namespace biofvm
} // namespace physicore
