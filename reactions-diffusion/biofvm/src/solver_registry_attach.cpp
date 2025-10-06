#include <cassert>

#ifdef PHYSICORE_HAS_THRUST
	#include "thrust_solver/register_solver.h"
#endif

#include "openmp_solver/register_solver.h"
#include "solver_registry_sole.cpp"

struct attachment_point
{
	attachment_point()
	{
		kernels::openmp_solver::attach_to_registry();
#ifdef PHYSICORE_HAS_THRUST
		kernels::thrust_solver::attach_to_registry();
#endif
	}
};

static const attachment_point ap;
