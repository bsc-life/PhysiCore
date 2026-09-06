#pragma once

#include <common/generic_agent_solver.h>
#include <physicell/environment.h>

namespace physicore::mechanics::physicell::kernels::openmp_solver {

class position_solver : private generic_agent_solver<mechanical_agent>
{
public:
	void update_cell_forces(environment& e);

	void update_cell_neighbors(environment& e, const cartesian_mesh& mesh);

	void update_motility(environment& e);

	void update_basement_membrane_interactions(environment& e, const cartesian_mesh& mesh);

	void update_spring_attachments(environment& e);

	void update_positions(environment& e);
};

} // namespace physicore::mechanics::physicell::kernels::openmp_solver
