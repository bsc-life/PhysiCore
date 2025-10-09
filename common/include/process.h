#pragma once

namespace physicore {

class timestep_executor
{
public:
	virtual void run_single_timestep() = 0;

	virtual void serialize_state() = 0;

	virtual ~timestep_executor() = default;
};

} // namespace physicore
