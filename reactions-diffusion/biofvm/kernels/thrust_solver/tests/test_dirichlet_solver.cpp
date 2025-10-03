#include <gtest/gtest.h>
#include <noarr/structures/interop/bag.hpp>

#include "diffusion_solver.h"
#include "dirichlet_solver.h"
#include "microenvironment.h"

using namespace physicore;
using namespace physicore::biofvm;

using namespace physicore::biofvm::kernels::thrust_solver;

static std::unique_ptr<microenvironment> default_microenv(cartesian_mesh mesh)
{
	real_t timestep = 5;
	index_t substrates_count = 2;

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

static void add_dirichlet_at(microenvironment& m, index_t substrates_count,
							 const std::vector<std::array<index_t, 3>>& indices, const std::vector<real_t>& values)
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

static void add_boundary_dirichlet(microenvironment& m, index_t substrates_count, index_t dim_idx, bool min,
								   real_t value)
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


TEST(ThrustDirichletSolverTest, Interior1D)
{
	cartesian_mesh mesh(1, { 0, 0, 0 }, { 100, 0, 0 }, { 20, 0, 0 });

	auto m = default_microenv(mesh);

	add_dirichlet_at(*m, m->substrates_count, { { 2, 0, 0 } }, { 41 });

	diffusion_solver d_s;
	dirichlet_solver s;
	data_manager mgr;

	d_s.initialize(*m, 1);
	s.initialize(*m);
	mgr.initialize(*m, d_s);
	s.solve(*m, d_s);
	mgr.transfer_to_host();

	auto dens_l = d_s.get_substrates_layout<1>();
	real_t* densities = mgr.substrate_densities;

	noarr::traverser(dens_l).for_dims<'x'>([&](auto t) {
		auto s = t.state();

		auto l = dens_l ^ noarr::fix(s);

		if (noarr::get_index<'x'>(s) == 2)
		{
			EXPECT_FLOAT_EQ(l | noarr::get_at<'s'>(densities, 0), 41);
		}
	});
}

TEST(ThrustDirichletSolverTest, Interior2D)
{
	cartesian_mesh mesh(2, { 0, 0, 0 }, { 60, 60, 0 }, { 20, 20, 0 });

	auto m = default_microenv(mesh);

	add_dirichlet_at(*m, m->substrates_count, { { 1, 1, 0 }, { 2, 0, 0 } }, { 10, 11 });

	diffusion_solver d_s;
	dirichlet_solver s;
	data_manager mgr;

	d_s.initialize(*m, 1);
	s.initialize(*m);
	mgr.initialize(*m, d_s);
	s.solve(*m, d_s);
	mgr.transfer_to_host();

	auto dens_l = d_s.get_substrates_layout<2>();
	real_t* densities = mgr.substrate_densities;

	// second substrate should not change
	noarr::traverser(dens_l).for_dims<'x', 'y'>([&](auto t) {
		auto s = t.state();

		auto l = dens_l ^ noarr::fix(s);

		auto [x, y] = noarr::get_indices<'x', 'y'>(s);

		if (x == 1 && y == 1)
		{
			EXPECT_FLOAT_EQ(l | noarr::get_at<'s'>(densities, 0), 10);
		}

		if (x == 2 && y == 0)
		{
			EXPECT_FLOAT_EQ(l | noarr::get_at<'s'>(densities, 0), 11);
		}
	});
}

TEST(ThrustDirichletSolverTest, Interior3D)
{
	cartesian_mesh mesh(3, { 0, 0, 0 }, { 60, 60, 60 }, { 20, 20, 20 });

	auto m = default_microenv(mesh);

	add_dirichlet_at(*m, m->substrates_count, { { 1, 1, 1 }, { 1, 0, 2 }, { 0, 1, 2 } }, { 1000, 1001, 1002 });

	diffusion_solver d_s;
	dirichlet_solver s;
	data_manager mgr;

	d_s.initialize(*m, 1);
	s.initialize(*m);
	mgr.initialize(*m, d_s);
	s.solve(*m, d_s);
	mgr.transfer_to_host();

	auto dens_l = d_s.get_substrates_layout<3>();
	real_t* densities = mgr.substrate_densities;

	// second substrate should not change
	noarr::traverser(dens_l).for_dims<'x', 'y', 'z'>([&](auto t) {
		auto s = t.state();

		auto l = dens_l ^ noarr::fix(s);

		auto [x, y, z] = noarr::get_indices<'x', 'y', 'z'>(s);

		if (x == 1 && y == 1 && z == 1)
		{
			EXPECT_FLOAT_EQ(l | noarr::get_at<'s'>(densities, 0), 1000);
		}

		if (x == 1 && y == 0 && z == 2)
		{
			EXPECT_FLOAT_EQ(l | noarr::get_at<'s'>(densities, 0), 1001);
		}

		if (x == 0 && y == 1 && z == 2)
		{
			EXPECT_FLOAT_EQ(l | noarr::get_at<'s'>(densities, 0), 1002);
		}
	});
}

TEST(ThrustDirichletSolverTest, Boundary1D)
{
	cartesian_mesh mesh(1, { 0, 0, 0 }, { 100, 100, 100 }, { 20, 20, 20 });

	auto m = default_microenv(mesh);

	add_boundary_dirichlet(*m, m->substrates_count, 0, true, 4);
	add_boundary_dirichlet(*m, m->substrates_count, 0, false, 5);

	diffusion_solver d_s;
	dirichlet_solver s;
	data_manager mgr;

	d_s.initialize(*m, 1);
	s.initialize(*m);
	mgr.initialize(*m, d_s);
	s.solve(*m, d_s);
	mgr.transfer_to_host();

	auto dens_l = d_s.get_substrates_layout<1>();
	auto densities = noarr::make_bag(dens_l, mgr.substrate_densities);

	for (index_t x = 0; x < m->mesh.grid_shape[0]; x++)
	{
		if (x == 0)
			EXPECT_FLOAT_EQ((densities.at<'x', 's'>(x, 0)), 4);
		else if (x == m->mesh.grid_shape[0] - 1)
			EXPECT_FLOAT_EQ((densities.at<'x', 's'>(x, 0)), 5);
		else
			EXPECT_FLOAT_EQ((densities.at<'x', 's'>(x, 0)), 1);

		EXPECT_FLOAT_EQ((densities.at<'x', 's'>(x, 1)), 1);
	}
}

TEST(ThrustDirichletSolverTest, Boundary2D)
{
	cartesian_mesh mesh(2, { 0, 0, 0 }, { 100, 100, 100 }, { 20, 20, 20 });

	auto m = default_microenv(mesh);

	add_boundary_dirichlet(*m, m->substrates_count, 0, true, 4);
	add_boundary_dirichlet(*m, m->substrates_count, 0, false, 5);
	add_boundary_dirichlet(*m, m->substrates_count, 1, true, 6);
	add_boundary_dirichlet(*m, m->substrates_count, 1, false, 7);

	diffusion_solver d_s;
	dirichlet_solver s;
	data_manager mgr;

	d_s.initialize(*m, 1);
	s.initialize(*m);
	mgr.initialize(*m, d_s);
	s.solve(*m, d_s);
	mgr.transfer_to_host();

	auto dens_l = d_s.get_substrates_layout<2>();
	auto densities = noarr::make_bag(dens_l, mgr.substrate_densities);

	for (index_t x = 0; x < m->mesh.grid_shape[0]; x++)
		for (index_t y = 0; y < m->mesh.grid_shape[1]; y++)
		{
			// y boundary overwrites x boundary
			if (x == 0 && y == 0)
				EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 6);
			else if (x == 0 && y == m->mesh.grid_shape[1] - 1)
				EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 7);
			else if (x == m->mesh.grid_shape[0] - 1 && y == 0)
				EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 6);
			else if (x == m->mesh.grid_shape[0] - 1 && y == m->mesh.grid_shape[1] - 1)
				EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 7);

			// x boundary
			else if (x == 0)
				EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 4);
			else if (x == m->mesh.grid_shape[0] - 1)
				EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 5);

			// y boundary
			else if (y == 0)
				EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 6);
			else if (y == m->mesh.grid_shape[1] - 1)
				EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 7);

			// interior
			else
				EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(x, y, 0)), 1);

			EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(x, y, 1)), 1);
		}
}

TEST(ThrustDirichletSolverTest, Boundary3D)
{
	cartesian_mesh mesh(3, { 0, 0, 0 }, { 100, 100, 100 }, { 20, 20, 20 });

	auto m = default_microenv(mesh);

	add_boundary_dirichlet(*m, m->substrates_count, 0, true, 4);
	add_boundary_dirichlet(*m, m->substrates_count, 0, false, 5);
	add_boundary_dirichlet(*m, m->substrates_count, 1, true, 6);
	add_boundary_dirichlet(*m, m->substrates_count, 1, false, 7);
	add_boundary_dirichlet(*m, m->substrates_count, 2, true, 8);
	add_boundary_dirichlet(*m, m->substrates_count, 2, false, 9);

	diffusion_solver d_s;
	dirichlet_solver s;
	data_manager mgr;

	d_s.initialize(*m, 1);
	s.initialize(*m);
	mgr.initialize(*m, d_s);
	s.solve(*m, d_s);
	mgr.transfer_to_host();

	auto dens_l = d_s.get_substrates_layout<3>();
	auto densities = noarr::make_bag(dens_l, mgr.substrate_densities);

	// with only interior z indices
	for (index_t x = 0; x < m->mesh.grid_shape[0]; x++)
		for (index_t y = 0; y < m->mesh.grid_shape[1]; y++)
			for (index_t z = 1; z < m->mesh.grid_shape[2] - 1; z++)
			{
				// y boundary overwrites x boundary
				if (x == 0 && y == 0)
					EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 6);
				else if (x == 0 && y == m->mesh.grid_shape[1] - 1)
					EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 7);
				else if (x == m->mesh.grid_shape[0] - 1 && y == 0)
					EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 6);
				else if (x == m->mesh.grid_shape[0] - 1 && y == m->mesh.grid_shape[1] - 1)
					EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 7);

				// x boundary
				else if (x == 0)
					EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 4);
				else if (x == m->mesh.grid_shape[0] - 1)
					EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 5);

				// y boundary
				else if (y == 0)
					EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 6);
				else if (y == m->mesh.grid_shape[1] - 1)
					EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 7);

				// interior
				else
					EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 0)), 1);

				EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, z, 1)), 1);
			}

	// with exterior z indices
	for (index_t x = 0; x < m->mesh.grid_shape[0]; x++)
		for (index_t y = 0; y < m->mesh.grid_shape[1]; y++)
		{
			EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, 0, 0)), 8);
			EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(x, y, m->mesh.grid_shape[2] - 1, 0)), 9);
		}
}
