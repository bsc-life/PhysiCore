#include <biofvm/microenvironment.h>
#include <gtest/gtest.h>
#include <noarr/structures/interop/bag.hpp>

#include "bulk_solver.h"
#include "diffusion_solver.h"

using namespace physicore;
using namespace physicore::reactions_diffusion::biofvm;

using namespace physicore::reactions_diffusion::biofvm::kernels::openmp_solver;

namespace {
std::unique_ptr<microenvironment> default_microenv(cartesian_mesh mesh)
{
	const real_t timestep = 0.01;
	const index_t substrates_count = 2;

	auto diff_coefs = std::make_unique<real_t[]>(2);
	diff_coefs[0] = 4;
	diff_coefs[1] = 2;
	auto decay_rates = std::make_unique<real_t[]>(2);
	decay_rates[0] = 5;
	decay_rates[1] = 3;

	auto initial_conds = std::make_unique<real_t[]>(2);
	initial_conds[0] = 10;
	initial_conds[1] = 1;

	auto m = std::make_unique<microenvironment>(mesh, substrates_count, timestep);
	m->diffusion_coefficients = std::move(diff_coefs);
	m->decay_rates = std::move(decay_rates);
	m->initial_conditions = std::move(initial_conds);

	return m;
}
} // namespace

struct test_functor : bulk_functor
{
	real_t supply_rates(index_t s, index_t x, index_t y, index_t z) override
	{
		if (x == 1 && y == 1 && z == 1 && s == 0)
			return 5;
		return 0;
	}

	real_t supply_target_densities(index_t s, index_t x, index_t y, index_t z) override
	{
		if (x == 1 && y == 1 && z == 1 && s == 0)
			return 6;
		return 0;
	}

	real_t uptake_rates(index_t s, index_t x, index_t y, index_t z) override
	{
		if (x == 1 && y == 1 && z == 1 && s == 0)
			return 7;
		return 0;
	}
};

TEST(BulkSolverTest, SimulateBulkSource)
{
	const cartesian_mesh mesh(3, { 0, 0, 0 }, { 100, 100, 100 }, { 20, 20, 20 });

	auto m = default_microenv(mesh);

	m->bulk_fnc = std::make_unique<test_functor>();

	diffusion_solver d_s;
	bulk_solver solver;

	d_s.prepare(*m, 1);
	d_s.initialize();
	solver.initialize(*m);

#pragma omp parallel
	solver.solve(*m, d_s);

	auto dens_l = d_s.get_substrates_layout<3>();
	real_t* densities = d_s.get_substrates_pointer();

	for (index_t x = 0; x < (dens_l | noarr::get_length<'x'>()); x++)
		for (index_t y = 0; y < (dens_l | noarr::get_length<'y'>()); y++)
			for (index_t z = 0; z < (dens_l | noarr::get_length<'z'>()); z++)
				for (index_t s = 0; s < (dens_l | noarr::get_length<'s'>()); s++)
				{
					auto idx = noarr::idx<'s', 'x', 'y', 'z'>(s, x, y, z);

					if (x == 1 && y == 1 && z == 1 && s == 0)
						EXPECT_NEAR(dens_l | noarr::get_at(densities, idx), 9.19643, 1e-4);
					else if (s == 1)
						EXPECT_DOUBLE_EQ(dens_l | noarr::get_at(densities, idx), 1);
					else
						EXPECT_DOUBLE_EQ(dens_l | noarr::get_at(densities, idx), 10);
				}
}
