#include "physicell/openmp_solver/register_solver.h"

#include <physicell/solver_registry.h>

#include "openmp_solver.h"

void physicore::mechanics::physicell::kernels::openmp_solver::attach_to_registry()
{
	static const physicore::mechanics::physicell::registry_adder<openmp_solver> openmp_solver_adder("openmp_solver");
}
