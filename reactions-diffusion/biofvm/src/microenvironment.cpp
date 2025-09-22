#include "microenvironment.h"

#include "agent_container.h"

using namespace physicore::biofvm;

microenvironment::microenvironment(cartesian_mesh mesh, index_t substrates_count, real_t timestep)
	: mesh(mesh), agents(mesh.dims, substrates_count), substrates_count(substrates_count), diffusion_timestep(timestep)
{}

void microenvironment::run_single_timestep()
{
	// TODO: Implement the logic for running a single timestep of the microenvironment simulation
}
