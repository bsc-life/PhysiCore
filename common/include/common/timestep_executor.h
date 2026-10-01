#pragma once

#include "types.h"

namespace physicore {

class timestep_executor
{
public:
	timestep_executor() = default;
	timestep_executor(const timestep_executor&) = delete;
	timestep_executor(timestep_executor&&) = delete;
	timestep_executor& operator=(const timestep_executor&) = delete;
	timestep_executor& operator=(timestep_executor&&) = delete;

	virtual void run_single_timestep() = 0;

	virtual void serialize_state(real_t current_time) = 0;

	virtual ~timestep_executor() = default;
};

} // namespace physicore
