#include "thrust_solver.h"

using namespace physicore;
using namespace physicore::biofvm::kernels::thrust_solver;

void thrust_solver::solve(biofvm::microenvironment& m, index_t iterations)
{
	if (!initialized)
	{
		d_solver.initialize(m, 1);

		dir_solver.initialize(m);

		b_solver.initialize(m);

		c_solver.initialize(m);

		mgr.initialize(m, d_solver);

		initialized = true;
	}

	for (index_t it = 0; it < iterations; it++)
	{
		d_solver.solve();

		dir_solver.solve(m, d_solver);

		b_solver.solve(m, d_solver);

		c_solver.simulate_secretion_and_uptake(m, d_solver, mgr, it == 0);
	}
}
