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

	position_solver mechanics_position_solver;

	for (index_t i = 0; i < iterations; ++i)
	{
		mechanics_position_solver.update_cell_neighbors(e, e.mesh);

		mechanics_position_solver.update_cell_forces(e);

		mechanics_position_solver.update_motility(e);

		mechanics_position_solver.update_basement_membrane_interactions(e, e.mesh);

		mechanics_position_solver.update_spring_attachments(e);

		mechanics_position_solver.update_positions(e);
	}
}

} // namespace physicore::mechanics::physicell::kernels::openmp_solver
