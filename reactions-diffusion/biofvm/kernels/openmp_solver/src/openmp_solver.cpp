#include "openmp_solver.h"

#include <biofvm/microenvironment.h>

#include "dirichlet_solver.h"

using namespace physicore;
using namespace physicore::reactions_diffusion::biofvm::kernels::openmp_solver;

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

		dirichlet_solver::solve(m, d_solver);

		b_solver.solve(m, d_solver);

		c_solver.simulate_secretion_and_uptake(m, d_solver, recompute_cells);
	}

	recompute_cells = false;
}

real_t openmp_solver::get_substrate_density(index_t s, index_t x, index_t y, index_t z) const
{
	auto dens_l = d_solver.get_substrates_layout<3>();
	const auto* densities = d_solver.get_substrates_pointer();

	return dens_l | noarr::get_at<'s', 'x', 'y', 'z'>(densities, s, x, y, z);
}

real_t& openmp_solver::get_substrate_density(index_t s, index_t x, index_t y, index_t z)
{
	auto dens_l = d_solver.get_substrates_layout<3>();
	auto* densities = d_solver.get_substrates_pointer();

	return dens_l | noarr::get_at<'s', 'x', 'y', 'z'>(densities, s, x, y, z);
}

std::array<real_t, 3> openmp_solver::get_substrate_gradient(const microenvironment& m, index_t s, index_t x, index_t y,
															index_t z) const
{
	std::array<real_t, 3> gradient = { 0.0, 0.0, 0.0 };

	std::array<index_t, 3> coords = { x, y, z };

	for (index_t dim = 0; dim < m.mesh.dims; ++dim)
	{
		auto lower = coords;
		auto upper = coords;

		real_t distance = 2.0 * (real_t)m.mesh.voxel_shape[dim];

		if (coords[dim] > 0)
			lower[dim] -= 1;
		else
			distance = (real_t)m.mesh.voxel_shape[dim];

		if (coords[dim] < m.mesh.grid_shape[dim] - 1)
			upper[dim] += 1;
		else
			distance = (real_t)m.mesh.voxel_shape[dim];

		gradient[dim] = (get_substrate_density(s, upper[0], upper[1], upper[2])
						 - get_substrate_density(s, lower[0], lower[1], lower[2]))
						/ distance;
	}

	return gradient;
}

void openmp_solver::reinitialize_dirichlet([[maybe_unused]] microenvironment& m)
{
	// OpenMP solver doesn't need to reinitialize Dirichlet conditions
	// since it accesses the microenvironment data directly
}

void openmp_solver::recompute_positional_data([[maybe_unused]] microenvironment& m) { recompute_cells = true; }
