#pragma once

#include <filesystem>
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

	static std::unique_ptr<environment> create_from_config(const std::filesystem::path& config_file);

	std::unique_ptr<mechanical_agent_container> agents;
	solver_ptr solver;
	serializer_ptr serializer;

	real_t mechanics_timestep;
	real_t simulation_time = 0.0;
	cartesian_mesh mesh;

	index_t agent_types_count;
	std::vector<std::string> agent_type_names;

	bool automated_spring_adhesion = true;
	bool virtual_wall_at_domain_edges = true;
};

} // namespace physicore::mechanics::physicell
