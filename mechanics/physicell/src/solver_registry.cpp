#include "physicell/solver_registry.h"

#include "physicell/openmp_solver/register_solver.h"

namespace physicore::mechanics::physicell {

solver_registry& solver_registry::instance()
{
	static solver_registry r;
	return r;
}

struct attachment_point
{
	attachment_point()
	{
		kernels::openmp_solver::attach_to_registry();

	}
};

static const attachment_point ap;

} // namespace physicore::mechanics::physicell
