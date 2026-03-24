#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>
#include <physicell/openmp_solver/position_solver.h>
#include <physicell/openmp_solver/register_solver.h>

#include "physicell/environment.h"

using namespace physicore::mechanics::physicell;
using physicore::index_t;
using physicore::real_t;

namespace {

constexpr real_t kTolerance = static_cast<real_t>(1e-6);

//Helpers for test orchestration
void add_agent(environment& env, std::initializer_list<real_t> pos, real_t radius = 1, index_t type = 0,
               std::uint8_t movable = 1)
{
	ASSERT_TRUE(env.agents != nullptr);
	env.agents->create();

	auto& data = env.get_agent_data();
	const index_t dims = data.base_data.dims;
	const index_t idx = data.agents_count - 1;

	std::vector<real_t> p(pos);
	p.resize(static_cast<std::size_t>(dims), 0);
	for (index_t d = 0; d < dims; ++d)
		data.base_data.positions[idx * dims + d] = p[d];

	data.radius[idx] = radius;
	data.state_data.agent_type_index[idx] = type;
	data.state_data.is_movable[idx] = movable;

	data.mechanics_data.cell_cell_repulsion_strength[idx] = 1;
	data.mechanics_data.cell_cell_adhesion_strength[idx] = 1;
	data.mechanics_data.relative_maximum_adhesion_distance[idx] = 1;

	for (index_t t = 0; t < std::max<index_t>(data.agent_types_count, 1); ++t)
		data.mechanics_data.cell_adhesion_affinities[idx * data.agent_types_count + t] = 1;
}

void clear_kinematics_and_pressure(environment& env)
{
	auto& data = env.get_agent_data();
	std::ranges::fill(data.velocity, static_cast<real_t>(0));
	std::ranges::fill(data.previous_velocity, static_cast<real_t>(0));
	std::ranges::fill(data.state_data.simple_pressure, static_cast<real_t>(0));
}

void connect_pair(environment& env)
{
	auto& data = env.get_agent_data();
	data.state_data.neighbors[0] = { 1 };
	data.state_data.neighbors[1] = { 0 };
}

void run_update_cell_forces(environment& env)
{
	clear_kinematics_and_pressure(env);
	kernels::openmp_solver::position_solver::update_cell_forces(env);
}

real_t velocity_component(environment& env, index_t agent, index_t dim)
{
	auto& data = env.get_agent_data();
	return data.velocity[agent * data.base_data.dims + dim];
}

void set_uniform_affinity(environment& env, real_t value)
{
	auto& data = env.get_agent_data();
	for (index_t agent = 0; agent < data.agents_count; ++agent)
		for (index_t type = 0; type < data.agent_types_count; ++type)
			data.mechanics_data.cell_adhesion_affinities[agent * data.agent_types_count + type] = value;
}

}  // namespace

TEST(SolvePair, RepulsiveForce1D_Overlapping)
{
	environment env(0.1, 1, 1, 1);
	add_agent(env, { 0 });
	add_agent(env, { 0.5 });

	auto& data = env.get_agent_data();
	data.mechanics_data.cell_cell_adhesion_strength[0] = 0;
	data.mechanics_data.cell_cell_adhesion_strength[1] = 0;
	connect_pair(env);

	run_update_cell_forces(env);

	EXPECT_LT(velocity_component(env, 0, 0), 0);
	EXPECT_GT(velocity_component(env, 1, 0), 0);
	EXPECT_NEAR(velocity_component(env, 0, 0) + velocity_component(env, 1, 0), 0, kTolerance);
	EXPECT_GT(data.state_data.simple_pressure[0], 0);
	EXPECT_GT(data.state_data.simple_pressure[1], 0);
}

TEST(SolvePair, NoForce1D_FarApart)
{
	environment env(0.1, 1, 1, 1);
	add_agent(env, { 0 });
	add_agent(env, { 5 });
	connect_pair(env);

	run_update_cell_forces(env);

	auto& data = env.get_agent_data();
	EXPECT_NEAR(velocity_component(env, 0, 0), 0, kTolerance);
	EXPECT_NEAR(velocity_component(env, 1, 0), 0, kTolerance);
	EXPECT_NEAR(data.state_data.simple_pressure[0], 0, kTolerance);
	EXPECT_NEAR(data.state_data.simple_pressure[1], 0, kTolerance);
}

TEST(SolvePair, AdhesiveForce1D_InAdhesionRange)
{
	environment env(0.1, 1, 1, 1);
	add_agent(env, { 0 });
	add_agent(env, { 3 });

	auto& data = env.get_agent_data();
	data.mechanics_data.relative_maximum_adhesion_distance[0] = 2;
	data.mechanics_data.relative_maximum_adhesion_distance[1] = 2;
	connect_pair(env);

	run_update_cell_forces(env);

	EXPECT_GT(velocity_component(env, 0, 0), 0);
	EXPECT_LT(velocity_component(env, 1, 0), 0);
	EXPECT_NEAR(velocity_component(env, 0, 0) + velocity_component(env, 1, 0), 0, kTolerance);
	EXPECT_NEAR(data.state_data.simple_pressure[0], 0, kTolerance);
	EXPECT_NEAR(data.state_data.simple_pressure[1], 0, kTolerance);
}

TEST(SolvePair, NewtonsThirdLaw1D_ForceSymmetry)
{
	environment env(0.1, 1, 1, 1);
	add_agent(env, { 0 });
	add_agent(env, { 0.5 });

	auto& data = env.get_agent_data();
	data.mechanics_data.cell_cell_adhesion_strength[0] = 0;
	data.mechanics_data.cell_cell_adhesion_strength[1] = 0;
	connect_pair(env);

	run_update_cell_forces(env);

	EXPECT_NEAR(velocity_component(env, 0, 0), -velocity_component(env, 1, 0), kTolerance);
}

TEST(SolvePair, SimplePressure1D_Accumulates)
{
	environment env(0.1, 1, 1, 1);
	add_agent(env, { 0 });
	add_agent(env, { 0.5 });

	auto& data = env.get_agent_data();
	data.mechanics_data.cell_cell_adhesion_strength[0] = 0;
	data.mechanics_data.cell_cell_adhesion_strength[1] = 0;
	connect_pair(env);

	run_update_cell_forces(env);

	EXPECT_GT(data.state_data.simple_pressure[0], 0);
	EXPECT_GT(data.state_data.simple_pressure[1], 0);
	EXPECT_NEAR(data.state_data.simple_pressure[0], data.state_data.simple_pressure[1], kTolerance);
}

TEST(SolvePair, ZeroRepulsion1D_NoRepulsiveForce)
{
	environment env(0.1, 1, 1, 1);
	add_agent(env, { 0 });
	add_agent(env, { 0.5 });

	auto& data = env.get_agent_data();
	data.mechanics_data.cell_cell_repulsion_strength[0] = 0;
	data.mechanics_data.cell_cell_repulsion_strength[1] = 0;
	data.mechanics_data.cell_cell_adhesion_strength[0] = 0;
	data.mechanics_data.cell_cell_adhesion_strength[1] = 0;
	connect_pair(env);

	run_update_cell_forces(env);

	EXPECT_NEAR(velocity_component(env, 0, 0), 0, kTolerance);
	EXPECT_NEAR(velocity_component(env, 1, 0), 0, kTolerance);
	EXPECT_GT(data.state_data.simple_pressure[0], 0);
	EXPECT_GT(data.state_data.simple_pressure[1], 0);
}

TEST(SolvePair, RepulsiveForce2D_Overlapping)
{
	environment env(0.1, 2, 1, 1);
	add_agent(env, { 0, 0 });
	add_agent(env, { 0.5, 0.5 });

	auto& data = env.get_agent_data();
	data.mechanics_data.cell_cell_adhesion_strength[0] = 0;
	data.mechanics_data.cell_cell_adhesion_strength[1] = 0;
	connect_pair(env);

	run_update_cell_forces(env);

	EXPECT_LT(velocity_component(env, 0, 0), 0);
	EXPECT_LT(velocity_component(env, 0, 1), 0);
	EXPECT_GT(velocity_component(env, 1, 0), 0);
	EXPECT_GT(velocity_component(env, 1, 1), 0);
	EXPECT_NEAR(velocity_component(env, 0, 0) + velocity_component(env, 1, 0), 0, kTolerance);
	EXPECT_NEAR(velocity_component(env, 0, 1) + velocity_component(env, 1, 1), 0, kTolerance);
}

TEST(SolvePair, AdhesiveForce2D_InAdhesionRange)
{
	environment env(0.1, 2, 1, 1);
	add_agent(env, { 0, 0 });
	add_agent(env, { 3, 4 });

	auto& data = env.get_agent_data();
	data.mechanics_data.relative_maximum_adhesion_distance[0] = 3;
	data.mechanics_data.relative_maximum_adhesion_distance[1] = 3;
	connect_pair(env);

	run_update_cell_forces(env);

	EXPECT_GT(velocity_component(env, 0, 0), 0);
	EXPECT_GT(velocity_component(env, 0, 1), 0);
	EXPECT_LT(velocity_component(env, 1, 0), 0);
	EXPECT_LT(velocity_component(env, 1, 1), 0);
}

TEST(SolvePair, NewtonsThirdLaw2D_ForceSymmetry)
{
	environment env(0.1, 2, 1, 1);
	add_agent(env, { 0, 0 });
	add_agent(env, { 0.5, 0.5 });

	auto& data = env.get_agent_data();
	data.mechanics_data.cell_cell_adhesion_strength[0] = 0;
	data.mechanics_data.cell_cell_adhesion_strength[1] = 0;
	connect_pair(env);

	run_update_cell_forces(env);

	EXPECT_NEAR(velocity_component(env, 0, 0), -velocity_component(env, 1, 0), kTolerance);
	EXPECT_NEAR(velocity_component(env, 0, 1), -velocity_component(env, 1, 1), kTolerance);
}

TEST(SolvePair, RepulsiveForce3D_Overlapping)
{
	environment env(0.1, 3, 1, 1);
	add_agent(env, { 0, 0, 0 });
	add_agent(env, { 0.5, 0.5, 0.5 });

	auto& data = env.get_agent_data();
	data.mechanics_data.cell_cell_adhesion_strength[0] = 0;
	data.mechanics_data.cell_cell_adhesion_strength[1] = 0;
	connect_pair(env);

	run_update_cell_forces(env);

	EXPECT_LT(velocity_component(env, 0, 0), 0);
	EXPECT_LT(velocity_component(env, 0, 1), 0);
	EXPECT_LT(velocity_component(env, 0, 2), 0);
	EXPECT_GT(velocity_component(env, 1, 0), 0);
	EXPECT_GT(velocity_component(env, 1, 1), 0);
	EXPECT_GT(velocity_component(env, 1, 2), 0);
}

TEST(SolvePair, AdhesiveForce3D_DirectionalAccuracy)
{
	environment env(0.1, 3, 1, 1);
	add_agent(env, { 0, 0, 0 });
	add_agent(env, { 1, 2, 2 });

	auto& data = env.get_agent_data();
	data.mechanics_data.relative_maximum_adhesion_distance[0] = 2;
	data.mechanics_data.relative_maximum_adhesion_distance[1] = 2;
	connect_pair(env);

	run_update_cell_forces(env);

	const real_t vx = velocity_component(env, 0, 0);
	const real_t vy = velocity_component(env, 0, 1);
	const real_t vz = velocity_component(env, 0, 2);
	EXPECT_GT(vx, 0);
	EXPECT_GT(vy, 0);
	EXPECT_GT(vz, 0);
	EXPECT_NEAR(vy / vx, static_cast<real_t>(2), kTolerance);
	EXPECT_NEAR(vz / vx, static_cast<real_t>(2), kTolerance);
}

TEST(SolvePair, NewtonsThirdLaw3D_ForceSymmetry)
{
	environment env(0.1, 3, 1, 1);
	add_agent(env, { 0, 0, 0 });
	add_agent(env, { 0.5, 0.5, 0.5 });

	auto& data = env.get_agent_data();
	data.mechanics_data.cell_cell_adhesion_strength[0] = 0;
	data.mechanics_data.cell_cell_adhesion_strength[1] = 0;
	connect_pair(env);

	run_update_cell_forces(env);

	EXPECT_NEAR(velocity_component(env, 0, 0), -velocity_component(env, 1, 0), kTolerance);
	EXPECT_NEAR(velocity_component(env, 0, 1), -velocity_component(env, 1, 1), kTolerance);
	EXPECT_NEAR(velocity_component(env, 0, 2), -velocity_component(env, 1, 2), kTolerance);
}

TEST(SolvePair, ZeroDistance1D_Minimum)
{
	environment env(0.1, 1, 1, 1);
	add_agent(env, { 0 });
	add_agent(env, { 0 });

	auto& data = env.get_agent_data();
	data.mechanics_data.cell_cell_adhesion_strength[0] = 0;
	data.mechanics_data.cell_cell_adhesion_strength[1] = 0;
	connect_pair(env);

	run_update_cell_forces(env);

	EXPECT_TRUE(std::isfinite(velocity_component(env, 0, 0)));
	EXPECT_TRUE(std::isfinite(velocity_component(env, 1, 0)));
	EXPECT_TRUE(std::isfinite(data.state_data.simple_pressure[0]));
	EXPECT_TRUE(std::isfinite(data.state_data.simple_pressure[1]));
}

TEST(SolvePair, DifferentRadii1D_AsymmetricForce)
{
	environment env(0.1, 1, 1, 1);
	add_agent(env, { 0 }, 1);
	add_agent(env, { 1.5 }, 2);

	auto& data = env.get_agent_data();
	data.mechanics_data.cell_cell_adhesion_strength[0] = 0;
	data.mechanics_data.cell_cell_adhesion_strength[1] = 0;
	connect_pair(env);

	run_update_cell_forces(env);

	EXPECT_LT(velocity_component(env, 0, 0), 0);
	EXPECT_GT(velocity_component(env, 1, 0), 0);
	EXPECT_NEAR(velocity_component(env, 0, 0) + velocity_component(env, 1, 0), 0, kTolerance);
	EXPECT_GT(data.state_data.simple_pressure[0], 0);
	EXPECT_GT(data.state_data.simple_pressure[1], 0);
}

TEST(SolvePair, ZeroAffinity1D_NoAdhesion)
{
	environment env(0.1, 1, 1, 1);
	add_agent(env, { 0 });
	add_agent(env, { 3 });

	auto& data = env.get_agent_data();
	data.mechanics_data.relative_maximum_adhesion_distance[0] = 2;
	data.mechanics_data.relative_maximum_adhesion_distance[1] = 2;
	set_uniform_affinity(env, 0);
	connect_pair(env);

	run_update_cell_forces(env);

	EXPECT_NEAR(velocity_component(env, 0, 0), 0, kTolerance);
	EXPECT_NEAR(velocity_component(env, 1, 0), 0, kTolerance);
	EXPECT_NEAR(data.state_data.simple_pressure[0], 0, kTolerance);
	EXPECT_NEAR(data.state_data.simple_pressure[1], 0, kTolerance);
}

TEST(SolvePair, DifferentCellTypes1D_AffinityLookup)
{
	environment env(0.1, 1, 2, 1);
	add_agent(env, { 0 }, 1, 0);
	add_agent(env, { 3 }, 1, 1);

	auto& data = env.get_agent_data();
	data.mechanics_data.relative_maximum_adhesion_distance[0] = 2;
	data.mechanics_data.relative_maximum_adhesion_distance[1] = 2;
	set_uniform_affinity(env, 0);
	data.mechanics_data.cell_adhesion_affinities[0 * data.agent_types_count + 1] = 1;
	data.mechanics_data.cell_adhesion_affinities[1 * data.agent_types_count + 0] = 1;
	connect_pair(env);

	run_update_cell_forces(env);

	EXPECT_GT(velocity_component(env, 0, 0), 0);
	EXPECT_LT(velocity_component(env, 1, 0), 0);
	EXPECT_NEAR(velocity_component(env, 0, 0) + velocity_component(env, 1, 0), 0, kTolerance);
}
