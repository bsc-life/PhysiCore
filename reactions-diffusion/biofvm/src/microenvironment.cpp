#include "microenvironment.h"

#include "agent_container.h"
#include "base_agent_data.h"

using namespace physicore::biofvm;

microenvironment::microenvironment(const cartesian_mesh& mesh, index_t substrates_count, real_t timestep)
	: diffusion_timestep(timestep), mesh(mesh), substrates_count(substrates_count)
{
	auto base_data = std::make_unique<base_agent_data>(mesh.dims);
	auto data = std::make_unique<agent_data>(*base_data, substrates_count);
	agents = make_unique<agent_container>(std::move(base_data), std::move(data));
}

void microenvironment::run_single_timestep() { solver->solve(*this, 1); }
