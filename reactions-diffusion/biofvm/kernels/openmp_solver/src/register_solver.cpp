#include "openmp_solver/register_solver.h"

#include "../../../include/solver_registry.h"
#include "openmp_solver.h"

void physicore::biofvm::kernels::openmp_solver::attach_to_registry()
{
	static physicore::biofvm::registry_adder<openmp_solver> openmp_solver_adder("openmp_solver");
}
