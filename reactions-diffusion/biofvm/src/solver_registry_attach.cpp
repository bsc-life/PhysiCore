#include <cassert>

#include "openmp_solver/register_solver.h"
#include "solver_registry_sole.cpp"
#include "thrust_solver/register_solver.h"

struct attachment_point
{
	attachment_point()
	{
		kernels::openmp_solver::attach_to_registry();
		kernels::thrust_solver::attach_to_registry();
	}
};

static const attachment_point ap;
