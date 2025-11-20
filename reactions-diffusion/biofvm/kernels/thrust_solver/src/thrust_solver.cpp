#include "thrust_solver.h"

using namespace physicore;
using namespace physicore::biofvm::kernels::thrust_solver;

void thrust_solver::initialize(biofvm::microenvironment& m)
{
	if (initialized)
		return;

	d_solver.initialize(m, 1);

	dir_solver.initialize(m);

	c_solver.initialize(m);

	mgr.initialize(m, d_solver);

	initialized = true;
}

void thrust_solver::solve(biofvm::microenvironment& m, index_t iterations)
{
	initialize(m);

	for (index_t it = 0; it < iterations; it++)
	{
		d_solver.solve();

		dir_solver.solve(m, d_solver);

		b_solver.solve(m, d_solver);

		c_solver.simulate_secretion_and_uptake(m, d_solver, mgr, it == 0);
	}
}

real_t thrust_solver::get_substrate_density(index_t s, index_t x, index_t y, index_t z) const
{
	auto dens_l = d_solver.get_substrates_layout<3>();
	auto densities = mgr.substrate_densities;

	return dens_l | noarr::get_at<'s', 'x', 'y', 'z'>(densities, s, x, y, z);
}
