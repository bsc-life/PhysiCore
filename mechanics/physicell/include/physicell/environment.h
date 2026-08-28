#pragma once

#include <memory>

#include <common/mesh.h>
#include <common/timestep_executor.h>
#include <common/types.h>

#include "mechanical_agent_container.h"
#include "serializer.h"
#include "solver.h"

namespace physicore::mechanics::physicell {

class environment : public timestep_executor
{
public:
	environment(const cartesian_mesh& mesh, index_t agent_types_count, index_t substrates_count, real_t timestep);

	void run_single_timestep() override;

	void serialize_state(real_t current_time) override;

	real_t mechanics_timestep;
	bool automated_spring_adhesion = true;
	bool virtual_wall_at_domain_edges = true;

	serializer_ptr serializer;
	solver_ptr solver;

	std::unique_ptr<mechanical_agent_container> agents;

	cartesian_mesh mesh;
};

} // namespace physicore::mechanics::physicell
