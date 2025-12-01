#include <biofvm/microenvironment.h>
#include <biofvm/microenvironment_builder.h>
#include <gtest/gtest.h>
#include <noarr/structures/interop/bag.hpp>

#include "diffusion_solver.h"
#include "dirichlet_solver.h"

using namespace physicore;
using namespace physicore::biofvm;

using namespace physicore::biofvm::kernels::openmp_solver;

namespace {
std::unique_ptr<microenvironment> default_microenv(cartesian_mesh mesh)
{
	const real_t timestep = 5;
	const index_t substrates_count = 2;

	auto diff_coefs = std::make_unique<real_t[]>(2);
	diff_coefs[0] = 4;
	diff_coefs[1] = 2;
	auto decay_rates = std::make_unique<real_t[]>(2);
	decay_rates[0] = 5;
	decay_rates[1] = 3;

	auto initial_conds = std::make_unique<real_t[]>(2);
	initial_conds[0] = 1;
	initial_conds[1] = 1;

	auto m = std::make_unique<microenvironment>(mesh, substrates_count, timestep);
	m->diffusion_coefficients = std::move(diff_coefs);
	m->decay_rates = std::move(decay_rates);
	m->initial_conditions = std::move(initial_conds);

	return m;
}

void add_dirichlet_at(microenvironment& m, index_t substrates_count, const std::vector<std::array<index_t, 3>>& indices,
					  const std::vector<real_t>& values)
{
	m.dirichlet_interior_voxels_count = indices.size();
	m.dirichlet_interior_voxels = std::make_unique<index_t[]>(m.dirichlet_interior_voxels_count * m.mesh.dims);

	for (index_t i = 0; i < m.dirichlet_interior_voxels_count; i++)
		for (index_t d = 0; d < m.mesh.dims; d++)
			m.dirichlet_interior_voxels[i * m.mesh.dims + d] = indices[i][d];

	m.dirichlet_interior_values = std::make_unique<real_t[]>(substrates_count * m.dirichlet_interior_voxels_count);
	m.dirichlet_interior_conditions = std::make_unique<bool[]>(substrates_count * m.dirichlet_interior_voxels_count);

	for (index_t i = 0; i < m.dirichlet_interior_voxels_count; i++)
	{
		m.dirichlet_interior_values[i * substrates_count] = values[i];
		m.dirichlet_interior_conditions[i * substrates_count] = true; // only the first substrate

		for (index_t j = 1; j < m.substrates_count; j++)
			m.dirichlet_interior_conditions[i * substrates_count + j] = false;
	}
}

void add_boundary_dirichlet(microenvironment& m, index_t substrates_count, index_t dim_idx, bool min, real_t value)
{
	auto& values = min ? m.dirichlet_min_boundary_values[dim_idx] : m.dirichlet_max_boundary_values[dim_idx];
	auto& conditions =
		min ? m.dirichlet_min_boundary_conditions[dim_idx] : m.dirichlet_max_boundary_conditions[dim_idx];

	if (!values)
	{
		values = std::make_unique<real_t[]>(substrates_count);
		conditions = std::make_unique<bool[]>(substrates_count);

		for (index_t s = 0; s < substrates_count; s++)
		{
			values[s] = 42;
			conditions[s] = false;
		}
	}

	// only the first substrate
	values[0] = value;
	conditions[0] = true;
}
} // namespace

TEST(DirichletSolverTest, Interior1D)
{
	const cartesian_mesh mesh(1, { 0, 0, 0 }, { 100, 0, 0 }, { 20, 0, 0 });

	auto m = default_microenv(mesh);

	add_dirichlet_at(*m, m->substrates_count, { { 2, 0, 0 } }, { 41 });

	diffusion_solver d_s;

	d_s.prepare(*m, 1);
	d_s.initialize();

#pragma omp parallel
	dirichlet_solver::solve(*m, d_s);

	auto dens_l = d_s.get_substrates_layout<1>();
	real_t* densities = d_s.get_substrates_pointer();

	noarr::traverser(dens_l).for_dims<'x'>([&](auto t) {
		auto s = t.state();

		auto l = dens_l ^ noarr::fix(s);

		if (noarr::get_index<'x'>(s) == 2)
		{
			EXPECT_DOUBLE_EQ(l | noarr::get_at<'s'>(densities, 0), 41);
		}
	});
}

TEST(DirichletSolverTest, Interior2D)
{
	const cartesian_mesh mesh(2, { 0, 0, 0 }, { 60, 60, 0 }, { 20, 20, 0 });

	auto m = default_microenv(mesh);

	add_dirichlet_at(*m, m->substrates_count, { { 1, 1, 0 }, { 2, 0, 0 } }, { 10, 11 });

	diffusion_solver d_s;

	d_s.prepare(*m, 1);
	d_s.initialize();

#pragma omp parallel
	dirichlet_solver::solve(*m, d_s);

	auto dens_l = d_s.get_substrates_layout<2>();
	real_t* densities = d_s.get_substrates_pointer();

	// second substrate should not change
	noarr::traverser(dens_l).for_dims<'x', 'y'>([&](auto t) {
		auto s = t.state();

		auto l = dens_l ^ noarr::fix(s);

		auto [x, y] = noarr::get_indices<'x', 'y'>(s);

		if (x == 1 && y == 1)
		{
			EXPECT_DOUBLE_EQ(l | noarr::get_at<'s'>(densities, 0), 10);
		}

		if (x == 2 && y == 0)
		{
			EXPECT_DOUBLE_EQ(l | noarr::get_at<'s'>(densities, 0), 11);
		}
	});
}

TEST(DirichletSolverTest, Interior3D)
{
	const cartesian_mesh mesh(3, { 0, 0, 0 }, { 60, 60, 60 }, { 20, 20, 20 });

	auto m = default_microenv(mesh);

	add_dirichlet_at(*m, m->substrates_count, { { 1, 1, 1 }, { 1, 0, 2 }, { 0, 1, 2 } }, { 1000, 1001, 1002 });

	diffusion_solver d_s;

	d_s.prepare(*m, 1);
	d_s.initialize();

#pragma omp parallel
	dirichlet_solver::solve(*m, d_s);

	auto dens_l = d_s.get_substrates_layout<3>();
	real_t* densities = d_s.get_substrates_pointer();

	// second substrate should not change
	noarr::traverser(dens_l).for_dims<'x', 'y', 'z'>([&](auto t) {
		auto s = t.state();

		auto l = dens_l ^ noarr::fix(s);

		auto [x, y, z] = noarr::get_indices<'x', 'y', 'z'>(s);

		if (x == 1 && y == 1 && z == 1)
		{
			EXPECT_DOUBLE_EQ(l | noarr::get_at<'s'>(densities, 0), 1000);
		}

		if (x == 1 && y == 0 && z == 2)
		{
			EXPECT_DOUBLE_EQ(l | noarr::get_at<'s'>(densities, 0), 1001);
		}

		if (x == 0 && y == 1 && z == 2)
		{
			EXPECT_DOUBLE_EQ(l | noarr::get_at<'s'>(densities, 0), 1002);
		}
	});
}

TEST(DirichletSolverTest, Boundary1D)
{
	const cartesian_mesh mesh(1, { 0, 0, 0 }, { 100, 100, 100 }, { 20, 20, 20 });

	auto m = default_microenv(mesh);

	add_boundary_dirichlet(*m, m->substrates_count, 0, true, 4);
	add_boundary_dirichlet(*m, m->substrates_count, 0, false, 5);

	diffusion_solver d_s;

	d_s.prepare(*m, 1);
	d_s.initialize();

#pragma omp parallel
	dirichlet_solver::solve(*m, d_s);

	auto dens_l = d_s.get_substrates_layout<1>();
	auto densities = noarr::make_bag(dens_l, d_s.get_substrates_pointer());

	for (index_t x = 0; x < m->mesh.grid_shape[0]; x++)
	{
		if (x == 0)
			EXPECT_DOUBLE_EQ((densities.at<'x', 's'>(x, 0)), 4);
		else if (x == m->mesh.grid_shape[0] - 1)
			EXPECT_DOUBLE_EQ((densities.at<'x', 's'>(x, 0)), 5);
		else
			EXPECT_DOUBLE_EQ((densities.at<'x', 's'>(x, 0)), 1);

		EXPECT_DOUBLE_EQ((densities.at<'x', 's'>(x, 1)), 1);
	}
}

TEST(DirichletSolverTest, Boundary2D)
{
	const cartesian_mesh mesh(2, { 0, 0, 0 }, { 100, 100, 100 }, { 20, 20, 20 });

	auto m = default_microenv(mesh);

	add_boundary_dirichlet(*m, m->substrates_count, 0, true, 4);
	add_boundary_dirichlet(*m, m->substrates_count, 0, false, 5);
	add_boundary_dirichlet(*m, m->substrates_count, 1, true, 6);
	add_boundary_dirichlet(*m, m->substrates_count, 1, false, 7);

	diffusion_solver d_s;

	d_s.prepare(*m, 1);
	d_s.initialize();

#pragma omp parallel
	dirichlet_solver::solve(*m, d_s);

	auto dens_l = d_s.get_substrates_layout<2>();
	auto densities = noarr::make_bag(dens_l, d_s.get_substrates_pointer());

	for (index_t x = 0; x < m->mesh.grid_shape[0]; x++)
		for (index_t y = 0; y < m->mesh.grid_shape[1]; y++)
		{
			// y boundary overwrites x boundary
			if (x == 0 && y == 0)
				EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 6);
			else if (x == 0 && y == m->mesh.grid_shape[1] - 1)
				EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 7);
			else if (x == m->mesh.grid_shape[0] - 1 && y == 0)
				EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 6);
			else if (x == m->mesh.grid_shape[0] - 1 && y == m->mesh.grid_shape[1] - 1)
				EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 7);

			// x boundary
			else if (x == 0)
				EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 4);
			else if (x == m->mesh.grid_shape[0] - 1)
				EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 5);

			// y boundary
			else if (y == 0)
				EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 6);
			else if (y == m->mesh.grid_shape[1] - 1)
				EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 7);

			// interior
			else
				EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 1);

			EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 's'>(x, y, 1)), 1);
		}
}

TEST(DirichletSolverTest, Boundary3D)
{
	const cartesian_mesh mesh(3, { 0, 0, 0 }, { 100, 100, 100 }, { 20, 20, 20 });

	auto m = default_microenv(mesh);

	add_boundary_dirichlet(*m, m->substrates_count, 0, true, 4);
	add_boundary_dirichlet(*m, m->substrates_count, 0, false, 5);
	add_boundary_dirichlet(*m, m->substrates_count, 1, true, 6);
	add_boundary_dirichlet(*m, m->substrates_count, 1, false, 7);
	add_boundary_dirichlet(*m, m->substrates_count, 2, true, 8);
	add_boundary_dirichlet(*m, m->substrates_count, 2, false, 9);

	diffusion_solver d_s;

	d_s.prepare(*m, 1);
	d_s.initialize();

#pragma omp parallel
	dirichlet_solver::solve(*m, d_s);

	auto dens_l = d_s.get_substrates_layout<3>();
	auto densities = noarr::make_bag(dens_l, d_s.get_substrates_pointer());

	// with only interior z indices
	for (index_t x = 0; x < m->mesh.grid_shape[0]; x++)
		for (index_t y = 0; y < m->mesh.grid_shape[1]; y++)
			for (index_t z = 1; z < m->mesh.grid_shape[2] - 1; z++)
			{
				// y boundary overwrites x boundary
				if (x == 0 && y == 0)
					EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 6);
				else if (x == 0 && y == m->mesh.grid_shape[1] - 1)
					EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 7);
				else if (x == m->mesh.grid_shape[0] - 1 && y == 0)
					EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 6);
				else if (x == m->mesh.grid_shape[0] - 1 && y == m->mesh.grid_shape[1] - 1)
					EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 7);

				// x boundary
				else if (x == 0)
					EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 4);
				else if (x == m->mesh.grid_shape[0] - 1)
					EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 5);

				// y boundary
				else if (y == 0)
					EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 6);
				else if (y == m->mesh.grid_shape[1] - 1)
					EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 7);

				// interior
				else
					EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 1);

				EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 1)), 1);
			}

	// with exterior z indices
	for (index_t x = 0; x < m->mesh.grid_shape[0]; x++)
		for (index_t y = 0; y < m->mesh.grid_shape[1]; y++)
		{
			EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, 0, 0)), 8);
			EXPECT_DOUBLE_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, m->mesh.grid_shape[2] - 1, 0)), 9);
		}
}

TEST(DirichletSolverTest, DynamicBoundaryModification)
{
	// Build a 3D microenvironment with oxygen substrate
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 100.0, 0.01, 38.0);
	builder.resize(3, { 0, 0, 0 }, { 100, 100, 100 }, { 10, 10, 10 });

	// Start with atmospheric oxygen at x-min boundary
	builder.add_boundary_dirichlet_conditions(0, // O2
											  { 160.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { true, false, false },
											  { false, false, false });

	auto env = builder.build();

	// Initialize the solver
	env->solver->initialize(*env);

	// Run a few timesteps
	for (int i = 0; i < 5; ++i)
	{
		env->run_single_timestep();
	}

	// Simulate a change in boundary conditions (e.g., hypoxic conditions)
	env->update_dirichlet_boundary_min('x', 0, 20.0, true); // Reduce oxygen at x-min boundary

	// Update the solver with new boundary conditions
	env->update_dirichlet_conditions();

	// Continue simulation with new boundary conditions
	for (int i = 0; i < 5; ++i)
	{
		env->run_single_timestep();
	}

	// Verify the boundary value was updated
	EXPECT_DOUBLE_EQ(env->get_substrate_density(0, 0, 1, 1), 20.0);
}

TEST(DirichletSolverTest, DynamicInteriorVoxelModification)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 100.0, 0.01, 38.0);
	builder.add_density("Glucose", "mM", 50.0, 0.02, 5.0);
	builder.resize(3, { 0, 0, 0 }, { 100, 100, 100 }, { 10, 10, 10 });

	// Add a constant source at the center
	builder.add_dirichlet_node({ 5, 5, 5 }, { 160.0, 50.0 }, { true, true });

	auto env = builder.build();
	env->solver->initialize(*env);

	// Run initial simulation
	for (int i = 0; i < 5; ++i)
	{
		env->run_single_timestep();
	}

	// Simulate depletion of the source (e.g., a nutrient reservoir running out)
	env->update_dirichlet_interior_voxel({ 5, 5, 5 }, 1, 0.0, true); // Deplete glucose at source

	env->update_dirichlet_conditions();

	// Continue simulation with depleted source
	for (int i = 0; i < 5; ++i)
	{
		env->run_single_timestep();
	}

	EXPECT_DOUBLE_EQ(env->get_substrate_density(1, 5, 5, 5), 0.0);
}

TEST(DirichletSolverTest, DisableDirichletConditionsDynamically)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 100.0, 0.01, 38.0);
	builder.resize(3, { 0, 0, 0 }, { 100, 100, 100 }, { 10, 10, 10 });

	// Start with boundary conditions
	builder.add_boundary_dirichlet_conditions(0, { 160.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { true, false, false },
											  { false, false, false });

	auto env = builder.build();
	env->solver->initialize(*env);

	// Run with boundary conditions
	for (int i = 0; i < 5; ++i)
	{
		env->run_single_timestep();
	}

	// Disable the boundary condition (e.g., removing a barrier)
	env->update_dirichlet_boundary_min('x', 0, 160.0, false);
	env->update_dirichlet_conditions();

	// Continue simulation without boundary condition
	for (int i = 0; i < 5; ++i)
	{
		env->run_single_timestep();
	}

	// Verify condition was disabled
	EXPECT_NE(env->get_substrate_density(1, 0, 1, 1), 160.0);
}

TEST(DirichletSolverTest, MultipleSubstrateModification)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 100.0, 0.01, 38.0);
	builder.add_density("Glucose", "mM", 50.0, 0.02, 5.0);
	builder.add_density("Lactate", "mM", 30.0, 0.015, 2.0);
	builder.resize(3, { 0, 0, 0 }, { 100, 100, 100 }, { 10, 10, 10 });

	// Add interior source with all three substrates
	builder.add_dirichlet_node({ 5, 5, 5 }, { 160.0, 50.0, 10.0 }, { true, true, true });

	auto env = builder.build();
	env->solver->initialize(*env);

	EXPECT_EQ(env->dirichlet_interior_voxels_count, 1);

	// Run initial simulation
	for (int i = 0; i < 5; ++i)
	{
		env->run_single_timestep();
	}

	env->update_dirichlet_interior_voxel({ 5, 5, 5 }, 0, 80.0, true);  // O2: reduce
	env->update_dirichlet_interior_voxel({ 5, 5, 5 }, 1, 100.0, true); // Glucose: increase
	env->update_dirichlet_interior_voxel({ 5, 5, 5 }, 2, 0.0, true);   // Lactate: deplete

	env->update_dirichlet_interior_voxel({ 2, 2, 2 }, 0, 200.0, true); // O2: increase
	env->update_dirichlet_interior_voxel({ 2, 2, 2 }, 1, 60.0, true);  // Glucose: increase
	env->update_dirichlet_interior_voxel({ 2, 2, 2 }, 2, 15.0, true);  // Lactate: increase

	env->update_dirichlet_conditions();

	// Continue simulation with modified voxels
	for (int i = 0; i < 5; ++i)
	{
		env->run_single_timestep();
	}

	EXPECT_DOUBLE_EQ(env->get_substrate_density(0, 5, 5, 5), 80.0);
	EXPECT_DOUBLE_EQ(env->get_substrate_density(1, 5, 5, 5), 100.0);
	EXPECT_DOUBLE_EQ(env->get_substrate_density(2, 5, 5, 5), 0.0);

	EXPECT_DOUBLE_EQ(env->get_substrate_density(0, 2, 2, 2), 200.0);
	EXPECT_DOUBLE_EQ(env->get_substrate_density(1, 2, 2, 2), 60.0);
	EXPECT_DOUBLE_EQ(env->get_substrate_density(2, 2, 2, 2), 15.0);
}
