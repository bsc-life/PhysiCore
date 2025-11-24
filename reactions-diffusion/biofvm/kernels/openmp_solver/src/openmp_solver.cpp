#include "openmp_solver.h"

using namespace physicore;
using namespace physicore::biofvm::kernels::openmp_solver;

void openmp_solver::initialize(biofvm::microenvironment& m)
{
	if (initialized)
		return;

	d_solver.prepare(m, 1);
	d_solver.initialize();

	b_solver.initialize(m);

	c_solver.initialize(m);

	initialized = true;
}

void openmp_solver::solve(biofvm::microenvironment& m, index_t iterations)
{
	initialize(m);

#pragma omp parallel
	for (index_t it = 0; it < iterations; it++)
	{
		d_solver.solve();

		dir_solver.solve(m, d_solver);

		b_solver.solve(m, d_solver);

		c_solver.simulate_secretion_and_uptake(m, d_solver, it == 0);
	}
}

real_t openmp_solver::get_substrate_density(index_t s, index_t x, index_t y, index_t z) const
{
	auto dens_l = d_solver.get_substrates_layout<3>();
	auto densities = d_solver.get_substrates_pointer();

	return dens_l | noarr::get_at<'s', 'x', 'y', 'z'>(densities, s, x, y, z);
}

real_t& openmp_solver::get_substrate_density(index_t s, index_t x, index_t y, index_t z)
{
	auto dens_l = d_solver.get_substrates_layout<3>();
	auto densities = d_solver.get_substrates_pointer();

	return dens_l | noarr::get_at<'s', 'x', 'y', 'z'>(densities, s, x, y, z);
}
