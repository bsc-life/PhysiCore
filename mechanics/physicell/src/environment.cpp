#include "environment.h"

#include <memory>

#include <common/base_agent_data.h>

using namespace physicore::mechanics::physicell;

environment::environment(const cartesian_mesh& mesh, index_t agent_types_count, index_t substrates_count,
						 real_t timestep)
	: mechanics_timestep(timestep), mesh(mesh)
{
	auto base_data = std::make_unique<physicore::base_agent_data>(mesh.dims);
	auto data = std::make_unique<mechanical_agent_data>(*base_data, agent_types_count, substrates_count);
	agents = std::make_unique<mechanical_agent_container>(std::move(base_data), std::move(data));
}

void environment::run_single_timestep() { solver->solve(*this, 1); }

void environment::serialize_state(real_t current_time)
{
	if (serializer)
	{
		serializer->serialize(current_time);
	}
}
