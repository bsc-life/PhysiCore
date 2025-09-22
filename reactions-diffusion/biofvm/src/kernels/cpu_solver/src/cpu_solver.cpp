#include "cpu_solver.h"

#include "diffusion_solver.h"

using namespace physicore;
using namespace physicore::biofvm::kernels::cpu;

void cpu_solver::solve(biofvm::microenvironment& m, index_t iterations)
{
	if (!initialized)
	{
		d_solver.prepare(m, 1);
		d_solver.initialize();

		b_solver.initialize(m);

		c_solver.initialize(m);

		initialized = true;
	}

#pragma omp parallel
	for (index_t it = 0; it < iterations; it++)
	{
		d_solver.solve();

		dir_solver.solve(m, d_solver);

		b_solver.solve(m, d_solver);

		c_solver.simulate_secretion_and_uptake(m, d_solver, it == 0);
	}
}
