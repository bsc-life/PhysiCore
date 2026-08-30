#pragma once

#include <optional>

#include <common/types.h>

#include "environment.h"

namespace physicore::mechanics::physicell {

class environment_builder
{
	real_t timestep = 0.1;
	real_t simulation_time = 0.0;
	std::optional<cartesian_mesh> mesh;

	std::vector<std::string> agent_type_names;

	std::vector<std::string> substrate_names;

	std::string solver_name = "openmp_solver";

public:
	void set_time_step(real_t time_step);
	void set_simulation_time(real_t sim_time);

	// mesh functions
	void resize(index_t dims, std::array<sindex_t, 3> bounding_box_mins, std::array<sindex_t, 3> bounding_box_maxs,
				std::array<index_t, 3> voxel_shape);

	void add_agent_type(std::string_view name);
	void add_substrate(std::string_view name);

	std::unique_ptr<environment> build();
};

} // namespace physicore::mechanics::physicell
