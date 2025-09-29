#include "microenvironment.h"

#include "agent_container.h"
#include "solver_registry.h"

using namespace physicore::biofvm;

microenvironment::microenvironment(const cartesian_mesh& mesh, index_t substrates_count, real_t timestep)
	: mesh(mesh),
	  agents(std::make_unique<agent_container>(mesh.dims, substrates_count)),
	  substrates_count(substrates_count),
	  diffusion_timestep(timestep)
{
	solver = solver_registry::instance().get("openmp_solver");
}

void microenvironment::run_single_timestep() { solver->solve(*this, 1); }
