#pragma once

#include <memory>

#include <common/timestep_executor.h>
#include <common/types.h>

#include "serializer.h"

namespace physicore::mechanics::physicell {

class environment : public timestep_executor
{
	real_t timestep;

public:
	explicit environment(real_t timestep) : timestep(timestep) {}

	void run_single_timestep() override;

	void serialize_state(real_t current_time) override;

	serializer_ptr serializer;
};

} // namespace physicore::mechanics::physicell
