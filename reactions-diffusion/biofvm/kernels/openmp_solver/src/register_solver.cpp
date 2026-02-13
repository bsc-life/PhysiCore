#include "openmp_solver/register_solver.h"

#include <biofvm/solver_registry.h>

#include "openmp_solver.h"

void physicore::reactions_diffusion::biofvm::kernels::openmp_solver::attach_to_registry()
{
	static const physicore::reactions_diffusion::biofvm::registry_adder<openmp_solver> openmp_solver_adder(
		"openmp_solver");
}
