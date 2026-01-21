#pragma once

#include <memory>

#include <common/cartesian_mesh.h>
#include <common/types.h>
#include <physicell/environment.h>


namespace physicore::mechanics::physicell::kernels::openmp_solver {

class position_solver
{
private:
public:
	// void prepare(const microenvironment& m, index_t iterations);

	// void initialize();

	// void solve();

	static void update_cell_forces(environment& e);

	static void update_cell_neighbors(environment& e, const cartesian_mesh& mesh);

	static void update_motility(environment& e);

	static void update_basement_membrane_interactions(environment& e, const cartesian_mesh& mesh);

	static void update_spring_attachments(environment& e);

	static void update_positions(environment& e);
};

} // namespace physicore::mechanics::physicell::kernels::openmp_solver
