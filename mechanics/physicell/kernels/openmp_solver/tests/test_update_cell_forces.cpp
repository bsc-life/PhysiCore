#include <algorithm>
#include <cmath>
#include <vector>

#include <gtest/gtest.h>
#include <physicell/openmp_solver/position_solver.h>
#include <physicell/openmp_solver/register_solver.h>

#include "physicell/environment.h"

using namespace physicore::mechanics::physicell;
using physicore::index_t;
using physicore::real_t;

namespace {

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
	std::fill(data.velocity.begin(), data.velocity.end(), static_cast<real_t>(0));
	std::fill(data.previous_velocity.begin(), data.previous_velocity.end(), static_cast<real_t>(0));
	std::fill(data.state_data.simple_pressure.begin(), data.state_data.simple_pressure.end(), static_cast<real_t>(0));
}

} // namespace

TEST(UpdateCellForcesTest, NoAgentsDoesNotCrash)
{
	environment env(0.1, 2, 1, 1);
	clear_kinematics_and_pressure(env);
	kernels::openmp_solver::position_solver::update_cell_forces(env);

	auto& data = env.get_agent_data();
	EXPECT_EQ(data.agents_count, 0);
	EXPECT_TRUE(data.velocity.empty());
	EXPECT_TRUE(data.state_data.simple_pressure.empty());
}

TEST(UpdateCellForcesTest, SingleAgentNoNeighborsLeavesZeroForces)
{
	environment env(0.1, 2, 1, 1);
	add_agent(env, { 0, 0 });
	clear_kinematics_and_pressure(env);

	kernels::openmp_solver::position_solver::update_cell_forces(env);

	auto& data = env.get_agent_data();
	EXPECT_EQ(data.agents_count, 1);
	EXPECT_FLOAT_EQ(data.velocity[0], 0);
	EXPECT_FLOAT_EQ(data.velocity[1], 0);
	EXPECT_FLOAT_EQ(data.state_data.simple_pressure[0], 0);
}

TEST(UpdateCellForcesTest, TwoAgentsRepelSymmetrically)
{
	environment env(0.1, 2, 1, 1);
	add_agent(env, { 0, 0 });
	add_agent(env, { 0.5, 0 });

	auto& data = env.get_agent_data();
	data.radius[0] = 1;
	data.radius[1] = 1;
	data.mechanics_data.cell_cell_adhesion_strength[0] = 0;
	data.mechanics_data.cell_cell_adhesion_strength[1] = 0;
	data.state_data.neighbors[0] = { 1 };
	data.state_data.neighbors[1] = { 0 };

	clear_kinematics_and_pressure(env);
	kernels::openmp_solver::position_solver::update_cell_forces(env);

	EXPECT_LT(data.velocity[0], 0);
	EXPECT_GT(data.velocity[2], 0);
	EXPECT_NEAR(data.velocity[0] + data.velocity[2], 0, 1e-6);
	EXPECT_NEAR(data.velocity[1], 0, 1e-6);
	EXPECT_NEAR(data.velocity[3], 0, 1e-6);
	EXPECT_GT(data.state_data.simple_pressure[0], 0);
	EXPECT_GT(data.state_data.simple_pressure[1], 0);
}

TEST(UpdateCellForcesTest, OverlappingAgentsProduceFiniteVelocities)
{
	environment env(0.1, 2, 1, 1);
	add_agent(env, { 0, 0 });
	add_agent(env, { 0, 0 });

	auto& data = env.get_agent_data();
	data.radius[0] = 1;
	data.radius[1] = 1;
	data.mechanics_data.cell_cell_adhesion_strength[0] = 0;
	data.mechanics_data.cell_cell_adhesion_strength[1] = 0;
	data.state_data.neighbors[0] = { 1 };

	clear_kinematics_and_pressure(env);
	kernels::openmp_solver::position_solver::update_cell_forces(env);

	for (real_t v : data.velocity)
		EXPECT_TRUE(std::isfinite(v));
	for (real_t p : data.state_data.simple_pressure)
		EXPECT_TRUE(std::isfinite(p));
}

TEST(UpdateCellForcesTest, AdhesionDependsOnAffinities)
{
	environment env(0.1, 2, 2, 1);
	add_agent(env, { 0, 0 }, 1, 0);
	add_agent(env, { 1, 0 }, 1, 1);

	auto& data = env.get_agent_data();
	data.mechanics_data.cell_cell_repulsion_strength[0] = 0;
	data.mechanics_data.cell_cell_repulsion_strength[1] = 0;
	data.mechanics_data.relative_maximum_adhesion_distance[0] = 2;
	data.mechanics_data.relative_maximum_adhesion_distance[1] = 2;
	data.state_data.neighbors[0] = { 1 };
	data.state_data.neighbors[1] = { 0 };

	// Both affinities enabled => attraction (agents move toward each other)
	data.mechanics_data.cell_adhesion_affinities[0 * 2 + 1] = 1;
	data.mechanics_data.cell_adhesion_affinities[1 * 2 + 0] = 1;
	clear_kinematics_and_pressure(env);
	kernels::openmp_solver::position_solver::update_cell_forces(env);
	EXPECT_GT(data.velocity[0], 0);
	EXPECT_LT(data.velocity[2], 0);

	// Disable one side affinity => no adhesion => no force
	data.mechanics_data.cell_adhesion_affinities[0 * 2 + 1] = 0;
	clear_kinematics_and_pressure(env);
	kernels::openmp_solver::position_solver::update_cell_forces(env);
	EXPECT_NEAR(data.velocity[0], 0, 1e-6);
	EXPECT_NEAR(data.velocity[2], 0, 1e-6);
}
