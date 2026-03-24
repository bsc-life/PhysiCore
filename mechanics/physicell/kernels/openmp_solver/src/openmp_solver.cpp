#include "openmp_solver.h"
#include <physicell/openmp_solver/position_solver.h>

namespace physicore::mechanics::physicell::kernels::openmp_solver {

void openmp_solver::initialize(environment& e)
{
	(void)e;
	initialized = true;
}

void openmp_solver::solve(environment& e, index_t iterations)
{

	if (!initialized)
	{
		initialize(e);
	}

	for (index_t i = 0; i < iterations; ++i)
	{
		position_solver::update_cell_neighbors(e, e.get_mesh());

		position_solver::update_cell_forces(e);

		position_solver::update_motility(e);

		position_solver::update_basement_membrane_interactions(e, e.get_mesh());

		position_solver::update_spring_attachments(e);

		position_solver::update_positions(e);
	}

}

} // namespace physicore::mechanics::physicell::kernels::openmp_solver
