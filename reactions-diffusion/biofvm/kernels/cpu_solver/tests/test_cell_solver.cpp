#include <gtest/gtest.h>
#include <noarr/structures/interop/bag.hpp>

#include "agent_data.h"
#include "cell_solver.h"
#include "diffusion_solver.h"
#include "generic_agent_solver.h"
#include "microenvironment.h"

using namespace physicore;
using namespace physicore::biofvm;

using namespace physicore::biofvm::kernels::cpu;

static std::unique_ptr<microenvironment> default_microenv(cartesian_mesh mesh, bool compute_interalized)
{
	real_t timestep = 0.01;
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

	m->compute_internalized_substrates = compute_interalized;

	return m;
}


static void compute_expected_agent_internalized_1d(auto densities, microenvironment& m, agent_data& agent_data,
												   std::vector<real_t>& expected_internalized)
{
	for (index_t i = 0; i < agent_data.base_data.agents_count; i++)
	{
		for (index_t s = 0; s < m.substrates_count; s++)
		{
			auto num = agent_data.get<secretion_uptake_component>()->secretion_rates[i * m.substrates_count + s]
					   * agent_data.get<secretion_uptake_component>()->saturation_densities[i * m.substrates_count + s]
					   * m.diffusion_timestep * agent_data.volumes[i];

			auto denom = (agent_data.get<secretion_uptake_component>()->secretion_rates[i * m.substrates_count + s]
						  + agent_data.get<secretion_uptake_component>()->uptake_rates[i * m.substrates_count + s])
						 * m.diffusion_timestep * agent_data.volumes[i] / m.mesh.voxel_volume();

			auto factor = agent_data.get<net_export_component>()->net_export_rates[i * m.substrates_count + s]
						  * m.diffusion_timestep;

			auto mesh_idx =
				m.mesh.voxel_position(std::span<const real_t>(agent_data.base_data.positions.data() + i, 1));

			expected_internalized[i * m.substrates_count + s] -=
				(m.mesh.voxel_volume() * -denom * densities.template at<'x', 's'>(mesh_idx[0], s) + num) / (1 + denom)
				+ factor;
		}
	}
}

static std::vector<real_t> compute_expected_agent_densities_1d(auto densities, microenvironment& m,
															   agent_data& agent_data)
{
	std::vector<real_t> expected_densities(m.mesh.voxel_count() * m.substrates_count, 0);

	for (index_t s = 0; s < m.substrates_count; s++)
	{
		for (index_t x = 0; x < m.mesh.grid_shape[0]; x++)
		{
			real_t num = 0, denom = 0, factor = 0;

			for (index_t i = 0; i < agent_data.base_data.agents_count; i++)
			{
				if (agent_data.base_data.positions[i] / m.mesh.voxel_shape[0] == x)
				{
					num +=
						agent_data.get<secretion_uptake_component>()->secretion_rates[i * m.substrates_count + s]
						* agent_data.get<secretion_uptake_component>()->saturation_densities[i * m.substrates_count + s]
						* m.diffusion_timestep * agent_data.volumes[i] / m.mesh.voxel_volume();

					denom += (agent_data.get<secretion_uptake_component>()->secretion_rates[i * m.substrates_count + s]
							  + agent_data.get<secretion_uptake_component>()->uptake_rates[i * m.substrates_count + s])
							 * m.diffusion_timestep * agent_data.volumes[i] / m.mesh.voxel_volume();

					factor += agent_data.get<net_export_component>()->net_export_rates[i * m.substrates_count + s]
							  * m.diffusion_timestep / m.mesh.voxel_volume();
				}
			}
			expected_densities[x * m.substrates_count + s] =
				((densities.template at<'x', 's'>(x, s) + num) / (1 + denom) + factor);
		}
	}

	return expected_densities;
}

void set_default_agent_values(agent* a, index_t rates_offset, index_t volume, std::array<real_t, 3> position,
							  index_t dims)
{
	a->secretion_rates()[0] = rates_offset + 100;
	a->secretion_rates()[1] = 0;

	a->uptake_rates()[0] = rates_offset + 200;
	a->uptake_rates()[1] = 0;

	a->saturation_densities()[0] = rates_offset + 300;
	a->saturation_densities()[1] = 0;

	a->net_export_rates()[0] = rates_offset + 400;
	a->net_export_rates()[1] = 0;

	a->fraction_released_at_death()[0] = 1;
	a->fraction_released_at_death()[1] = 0;

	a->volume() = volume;

	for (index_t i = 0; i < dims; ++i)
		a->position()[i] = position[i];
}

class RecomputeTest : public testing::TestWithParam<std::tuple<bool, bool>>
{};

INSTANTIATE_TEST_SUITE_P(CellSolverTest, RecomputeTest,
						 testing::Combine(testing::Values(true, false), testing::Values(true, false)));

TEST_P(RecomputeTest, Simple1D)
{
	bool compute_internalized = std::get<0>(GetParam());
	bool recompute = std::get<1>(GetParam());

	cartesian_mesh mesh(1, { 0, 0, 0 }, { 60, 20, 20 }, { 20, 20, 20 });

	auto m = default_microenv(mesh, compute_internalized);

	auto a1 = m->agents.create();
	auto a2 = m->agents.create();
	auto a3 = m->agents.create();

	set_default_agent_values(a1, 0, 1000, { 10, 0, 0 }, 1);
	set_default_agent_values(a2, 400, 1000, { 30, 0, 0 }, 1);
	set_default_agent_values(a3, 800, 1000, { 50, 0, 0 }, 1);

	diffusion_solver d_s;
	cell_solver s;

	d_s.prepare(*m, 1);
	d_s.initialize();
	s.initialize(*m);

#pragma omp parallel
	s.simulate_secretion_and_uptake(*m, d_s, true);

	auto dens_l = d_s.get_substrates_layout<1>();
	auto densities = noarr::make_bag(dens_l, d_s.get_substrates_pointer());

	if (compute_internalized)
	{
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[0], -216004.000000);
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[1], 0);

		EXPECT_FLOAT_EQ(a2->internalized_substrates()[0], -1469060.631579);
		EXPECT_FLOAT_EQ(a2->internalized_substrates()[1], 0);

		EXPECT_FLOAT_EQ(a3->internalized_substrates()[0], -2927715.703704);
		EXPECT_FLOAT_EQ(a3->internalized_substrates()[1], 0);
	}

	EXPECT_FLOAT_EQ((densities.template at<'x', 's'>(0, 0)), 28.000500);
	EXPECT_FLOAT_EQ((densities.template at<'x', 's'>(0, 1)), 1);

	EXPECT_FLOAT_EQ((densities.template at<'x', 's'>(1, 0)), 184.632579);
	EXPECT_FLOAT_EQ((densities.template at<'x', 's'>(1, 1)), 1);

	EXPECT_FLOAT_EQ((densities.template at<'x', 's'>(2, 0)), 366.964463);
	EXPECT_FLOAT_EQ((densities.template at<'x', 's'>(2, 1)), 1);

#pragma omp parallel
	s.simulate_secretion_and_uptake(*m, d_s, recompute);

	if (compute_internalized)
	{
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[0], -216004.000000 + -157093.818182);
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[1], 0);

		EXPECT_FLOAT_EQ(a2->internalized_substrates()[0], -1469060.631579 + -618551.844632);
		EXPECT_FLOAT_EQ(a2->internalized_substrates()[1], 0);

		EXPECT_FLOAT_EQ(a3->internalized_substrates()[0], -2927715.703704 + -867471.319407);
		EXPECT_FLOAT_EQ(a3->internalized_substrates()[1], 0);
	}

	EXPECT_FLOAT_EQ((densities.at<'x', 's'>(0, 0)), 47.637227);
	EXPECT_FLOAT_EQ((densities.at<'x', 's'>(0, 1)), 1);

	EXPECT_FLOAT_EQ((densities.at<'x', 's'>(1, 0)), 261.951560);
	EXPECT_FLOAT_EQ((densities.at<'x', 's'>(1, 1)), 1);

	EXPECT_FLOAT_EQ((densities.at<'x', 's'>(2, 0)), 475.398378);
	EXPECT_FLOAT_EQ((densities.at<'x', 's'>(2, 1)), 1);

	s.release_internalized_substrates(*m, d_s, 0);

	if (compute_internalized)
	{
		EXPECT_FLOAT_EQ((densities.at<'x', 's'>(0, 0)), 1);
		EXPECT_FLOAT_EQ((densities.at<'x', 's'>(0, 1)), 1);
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[0], 0);
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[1], 0);
	}
}

TEST_P(RecomputeTest, Simple2D)
{
	bool compute_internalized = std::get<0>(GetParam());
	bool recompute = std::get<1>(GetParam());

	cartesian_mesh mesh(2, { 0, 0, 0 }, { 60, 60, 20 }, { 20, 20, 20 });

	auto m = default_microenv(mesh, compute_internalized);

	auto a1 = m->agents.create();
	auto a2 = m->agents.create();
	auto a3 = m->agents.create();

	set_default_agent_values(a1, 0, 1000, { 10, 10, 0 }, 2);
	set_default_agent_values(a2, 400, 1000, { 30, 30, 0 }, 2);
	set_default_agent_values(a3, 800, 1000, { 50, 50, 0 }, 2);

	diffusion_solver d_s;
	cell_solver s;

	d_s.prepare(*m, 1);
	d_s.initialize();
	s.initialize(*m);

#pragma omp parallel
	s.simulate_secretion_and_uptake(*m, d_s, true);

	auto dens_l = d_s.get_substrates_layout<2>();
	auto densities = noarr::make_bag(dens_l, d_s.get_substrates_pointer());

	if (compute_internalized)
	{
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[0], -216004.000000);
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[1], 0);

		EXPECT_FLOAT_EQ(a2->internalized_substrates()[0], -1469060.631579);
		EXPECT_FLOAT_EQ(a2->internalized_substrates()[1], 0);

		EXPECT_FLOAT_EQ(a3->internalized_substrates()[0], -2927715.703704);
		EXPECT_FLOAT_EQ(a3->internalized_substrates()[1], 0);
	}

	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(0, 0, 0)), 28.000500);
	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(0, 0, 1)), 1);

	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(1, 1, 0)), 184.632579);
	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(1, 1, 1)), 1);

	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(2, 2, 0)), 366.964463);
	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(2, 2, 1)), 1);

#pragma omp parallel
	s.simulate_secretion_and_uptake(*m, d_s, recompute);

	if (compute_internalized)
	{
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[0], -216004.000000 + -157093.818182);
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[1], 0);

		EXPECT_FLOAT_EQ(a2->internalized_substrates()[0], -1469060.631579 + -618551.844632);
		EXPECT_FLOAT_EQ(a2->internalized_substrates()[1], 0);

		EXPECT_FLOAT_EQ(a3->internalized_substrates()[0], -2927715.703704 + -867471.319407);
		EXPECT_FLOAT_EQ(a3->internalized_substrates()[1], 0);
	}

	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(0, 0, 0)), 47.637227);
	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(0, 0, 1)), 1);

	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(1, 1, 0)), 261.951560);
	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(1, 1, 1)), 1);

	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(2, 2, 0)), 475.398378);
	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(2, 2, 1)), 1);

	s.release_internalized_substrates(*m, d_s, 0);

	if (compute_internalized)
	{
		EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(0, 0, 0)), 1);
		EXPECT_FLOAT_EQ((densities.at<'x', 'y', 's'>(0, 0, 1)), 1);
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[0], 0);
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[1], 0);
	}
}

TEST_P(RecomputeTest, Simple3D)
{
	bool compute_internalized = std::get<0>(GetParam());
	bool recompute = std::get<1>(GetParam());

	cartesian_mesh mesh(3, { 0, 0, 0 }, { 60, 60, 60 }, { 20, 20, 20 });

	auto m = default_microenv(mesh, compute_internalized);

	auto a1 = m->agents.create();
	auto a2 = m->agents.create();
	auto a3 = m->agents.create();

	set_default_agent_values(a1, 0, 1000, { 10, 10, 10 }, 3);
	set_default_agent_values(a2, 400, 1000, { 30, 30, 30 }, 3);
	set_default_agent_values(a3, 800, 1000, { 50, 50, 50 }, 3);

	diffusion_solver d_s;
	cell_solver s;

	d_s.prepare(*m, 1);
	d_s.initialize();
	s.initialize(*m);

#pragma omp parallel
	s.simulate_secretion_and_uptake(*m, d_s, true);

	auto dens_l = d_s.get_substrates_layout<3>();
	auto densities = noarr::make_bag(dens_l, d_s.get_substrates_pointer());

	if (compute_internalized)
	{
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[0], -216004.000000);
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[1], 0);

		EXPECT_FLOAT_EQ(a2->internalized_substrates()[0], -1469060.631579);
		EXPECT_FLOAT_EQ(a2->internalized_substrates()[1], 0);

		EXPECT_FLOAT_EQ(a3->internalized_substrates()[0], -2927715.703704);
		EXPECT_FLOAT_EQ(a3->internalized_substrates()[1], 0);
	}

	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(0, 0, 0, 0)), 28.000500);
	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(0, 0, 0, 1)), 1);

	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(1, 1, 1, 0)), 184.632579);
	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(1, 1, 1, 1)), 1);

	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(2, 2, 2, 0)), 366.964463);
	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(2, 2, 2, 1)), 1);

#pragma omp parallel
	s.simulate_secretion_and_uptake(*m, d_s, recompute);

	if (compute_internalized)
	{
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[0], -216004.000000 + -157093.818182);
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[1], 0);

		EXPECT_FLOAT_EQ(a2->internalized_substrates()[0], -1469060.631579 + -618551.844632);
		EXPECT_FLOAT_EQ(a2->internalized_substrates()[1], 0);

		EXPECT_FLOAT_EQ(a3->internalized_substrates()[0], -2927715.703704 + -867471.319407);
		EXPECT_FLOAT_EQ(a3->internalized_substrates()[1], 0);
	}

	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(0, 0, 0, 0)), 47.637227);
	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(0, 0, 0, 1)), 1);

	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(1, 1, 1, 0)), 261.951560);
	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(1, 1, 1, 1)), 1);

	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(2, 2, 2, 0)), 475.398378);
	EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(2, 2, 2, 1)), 1);

	s.release_internalized_substrates(*m, d_s, 0);

	if (compute_internalized)
	{
		EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(0, 0, 0, 0)), 1);
		EXPECT_FLOAT_EQ((densities.at<'x', 'y', 'z', 's'>(0, 0, 0, 1)), 1);
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[0], 0);
		EXPECT_FLOAT_EQ(a1->internalized_substrates()[1], 0);
	}
}

class agent_retriever : public generic_agent_solver<agent_data>
{};

TEST_P(RecomputeTest, Conflict)
{
	bool compute_internalized = std::get<0>(GetParam());
	bool recompute = std::get<1>(GetParam());

	cartesian_mesh mesh(1, { 0, 0, 0 }, { 60, 20, 20 }, { 20, 20, 20 });

	auto m = default_microenv(mesh, compute_internalized);

	std::vector<agent*> agents;

	for (int i = 0; i < 6; i++)
		agents.push_back(m->agents.create());

	set_default_agent_values(agents[0], 0, 500, { 10, 0, 0 }, 1);

	set_default_agent_values(agents[1], 600, 1000, { 30, 0, 0 }, 1);
	set_default_agent_values(agents[2], 1100, 1500, { 30, 0, 0 }, 1);

	set_default_agent_values(agents[3], 1600, 2000, { 50, 0, 0 }, 1);
	set_default_agent_values(agents[4], 2100, 2500, { 50, 0, 0 }, 1);
	set_default_agent_values(agents[5], 2600, 3000, { 50, 0, 0 }, 1);

	diffusion_solver d_s;
	cell_solver s;

	d_s.prepare(*m, 1);
	d_s.initialize();
	s.initialize(*m);

	auto dens_l = d_s.get_substrates_layout<1>();
	auto densities = noarr::make_bag(dens_l, d_s.get_substrates_pointer());

	auto& agent_data = agent_retriever().retrieve_agent_data(m->agents);

	std::vector<real_t> expected_internalized(agent_data.base_data.agents_count * m->substrates_count, 0);

	compute_expected_agent_internalized_1d(densities, *m, agent_data, expected_internalized);

#pragma omp parallel
	s.simulate_secretion_and_uptake(*m, d_s, true);

	if (compute_internalized)
	{
		for (std::size_t i = 0; i < agents.size(); i++)
		{
			EXPECT_FLOAT_EQ(agents[i]->internalized_substrates()[0], expected_internalized[2 * i]);
			EXPECT_FLOAT_EQ(agents[i]->internalized_substrates()[1], expected_internalized[2 * i + 1]);
		}
	}

	{
		auto expected = compute_expected_agent_densities_1d(densities, *m, agent_data);

		for (index_t x = 0; x < m->mesh.grid_shape[0]; x++)
		{
			EXPECT_FLOAT_EQ((densities.at<'x', 's'>(x, 0)), expected[2 * x]);
			EXPECT_FLOAT_EQ((densities.at<'x', 's'>(x, 1)), expected[2 * x + 1]);
		}
	}

	compute_expected_agent_internalized_1d(densities, *m, agent_data, expected_internalized);

#pragma omp parallel
	s.simulate_secretion_and_uptake(*m, d_s, recompute);

	if (compute_internalized)
	{
		for (std::size_t i = 0; i < agents.size(); i++)
		{
			EXPECT_FLOAT_EQ(agents[i]->internalized_substrates()[0], expected_internalized[2 * i]);
			EXPECT_FLOAT_EQ(agents[i]->internalized_substrates()[1], expected_internalized[2 * i + 1]);
		}
	}

	{
		auto expected = compute_expected_agent_densities_1d(densities, *m, agent_data);

		for (index_t x = 0; x < m->mesh.grid_shape[0]; x++)
		{
			EXPECT_FLOAT_EQ((densities.at<'x', 's'>(x, 0)), expected[2 * x]);
			EXPECT_FLOAT_EQ((densities.at<'x', 's'>(x, 1)), expected[2 * x + 1]);
		}
	}

#pragma omp parallel for
	for (std::size_t i = 0; i < agents.size(); i++)
	{
		s.release_internalized_substrates(*m, d_s, i);
	}

	if (compute_internalized)
	{
		for (index_t x = 0; x < m->mesh.grid_shape[0]; x++)
		{
			EXPECT_FLOAT_EQ((densities.at<'x', 's'>(x, 0)), 1);
			EXPECT_FLOAT_EQ((densities.at<'x', 's'>(x, 1)), 1);
		}
	}
}

TEST_P(RecomputeTest, ConflictBig)
{
	bool compute_internalized = std::get<0>(GetParam());
	bool recompute = std::get<1>(GetParam());
	index_t conflict_in_each_voxel = 50;

	cartesian_mesh mesh(1, { 0, 0, 0 }, { 2000, 20, 20 }, { 20, 20, 20 });

	auto m = default_microenv(mesh, compute_internalized);

	std::vector<agent*> agents;

	for (index_t i = 0; i < mesh.grid_shape[0]; i++)
	{
		for (index_t j = 0; j < conflict_in_each_voxel; j++)
		{
			agents.push_back(m->agents.create());
			set_default_agent_values(agents.back(), 0, 500, mesh.voxel_center({ i, 0, 0 }), 1);
		}
	}

	diffusion_solver d_s;
	cell_solver s;

	d_s.prepare(*m, 1);
	d_s.initialize();
	s.initialize(*m);

	auto dens_l = d_s.get_substrates_layout<1>();
	auto densities = noarr::make_bag(dens_l, d_s.get_substrates_pointer());

	auto& agent_data = agent_retriever().retrieve_agent_data(m->agents);

	std::vector<real_t> expected_internalized(agent_data.base_data.agents_count * m->substrates_count, 0);

	compute_expected_agent_internalized_1d(densities, *m, agent_data, expected_internalized);

#pragma omp parallel
	s.simulate_secretion_and_uptake(*m, d_s, true);

	if (compute_internalized)
	{
		for (std::size_t i = 0; i < agents.size(); i++)
		{
			EXPECT_FLOAT_EQ(agents[i]->internalized_substrates()[0], expected_internalized[2 * i]);
			EXPECT_FLOAT_EQ(agents[i]->internalized_substrates()[1], expected_internalized[2 * i + 1]);
		}
	}

	{
		auto expected = compute_expected_agent_densities_1d(densities, *m, agent_data);

		for (index_t x = 0; x < m->mesh.grid_shape[0]; x++)
		{
			EXPECT_FLOAT_EQ((densities.at<'x', 's'>(x, 0)), expected[2 * x]);
			EXPECT_FLOAT_EQ((densities.at<'x', 's'>(x, 1)), expected[2 * x + 1]);
		}
	}

	compute_expected_agent_internalized_1d(densities, *m, agent_data, expected_internalized);

#pragma omp parallel
	s.simulate_secretion_and_uptake(*m, d_s, recompute);

	if (compute_internalized)
	{
		for (std::size_t i = 0; i < agents.size(); i++)
		{
			EXPECT_FLOAT_EQ(agents[i]->internalized_substrates()[0], expected_internalized[2 * i]);
			EXPECT_FLOAT_EQ(agents[i]->internalized_substrates()[1], expected_internalized[2 * i + 1]);
		}
	}

	{
		auto expected = compute_expected_agent_densities_1d(densities, *m, agent_data);

		for (index_t x = 0; x < m->mesh.grid_shape[0]; x++)
		{
			EXPECT_FLOAT_EQ((densities.at<'x', 's'>(x, 0)), expected[2 * x]);
			EXPECT_FLOAT_EQ((densities.at<'x', 's'>(x, 1)), expected[2 * x + 1]);
		}
	}
}
