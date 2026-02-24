#pragma once

#include <memory>

#include <common/timestep_executor.h>
#include <common/types.h>

#include "mechanical_agent_container.h"
#include "serializer.h"
#include "solver.h"

namespace physicore::mechanics::physicell {

class environment : public timestep_executor
{
public:
	explicit environment(real_t timestep, index_t agent_types_count = 1, index_t substrates_count = 1)
		: environment(timestep, 3, agent_types_count, substrates_count)
	{}

	environment(real_t timestep, index_t dims, index_t agent_types_count, index_t substrates_count);

	void run_single_timestep() override;

	void serialize_state(real_t current_time) override;

	mechanical_agent_data& get_agent_data();

	real_t timestep;
	bool automated_spring_adhesion;
	bool virtual_wall_at_domain_edges;

	serializer_ptr serializer;
	solver_ptr solver;

	std::unique_ptr<mechanical_agent_container> agents;
};

} // namespace physicore::mechanics::physicell
