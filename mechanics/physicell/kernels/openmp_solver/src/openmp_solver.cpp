#include "openmp_solver.h"

namespace physicore::mechanics::physicell::kernels::openmp_solver {

void openmp_solver::initialize(environment& e)
{
	(void)e;
	initialized = true;
}

void openmp_solver::solve(environment& e, index_t iterations)
{
	(void)e;
	(void)iterations;
	if (!initialized)
	{
		initialize(e);
	}
}

} // namespace physicore::mechanics::physicell::kernels::openmp_solver
