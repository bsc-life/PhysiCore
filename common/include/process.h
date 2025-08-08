#pragma once

namespace physicore {

class process
{
public:
	virtual void run_single_timestep() = 0;

	virtual ~process() = default;
};

} // namespace physicore
