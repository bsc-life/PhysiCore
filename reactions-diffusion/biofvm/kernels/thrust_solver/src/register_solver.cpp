#include "thrust_solver/register_solver.h"

#include <biofvm/solver_registry.h>

#include "thrust_solver.h"

void physicore::biofvm::kernels::thrust_solver::attach_to_registry()
{
	static physicore::biofvm::registry_adder<thrust_solver> openmp_solver_adder("thrust_solver");
}
