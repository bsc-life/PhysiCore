#include <algorithm>
#include <vector>

#include <gtest/gtest.h>
#include <physicell/openmp_solver/position_solver.h>
#include <physicell/openmp_solver/register_solver.h>

#include "physicell/environment.h"

using namespace physicore::mechanics::physicell;
using physicore::cartesian_mesh;
using physicore::index_t;
using physicore::real_t;
using physicore::sindex_t;

namespace {

void add_agent(environment& env, std::initializer_list<real_t> pos, real_t radius = 1, std::uint8_t movable = 1,
			   real_t rel_max_adhesion_dist = 1)
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
	data.state_data.is_movable[idx] = movable;
	data.mechanics_data.relative_maximum_adhesion_distance[idx] = rel_max_adhesion_dist;

	data.mechanics_data.cell_cell_repulsion_strength[idx] = 1;
	data.mechanics_data.cell_cell_adhesion_strength[idx] = 1;

	for (index_t t = 0; t < std::max<index_t>(data.agent_types_count, 1); ++t)
		data.mechanics_data.cell_adhesion_affinities[idx * data.agent_types_count + t] = 1;
}

std::vector<index_t> sorted_neighbors(environment& env, index_t i)
{
	auto neighbors = env.get_agent_data().state_data.neighbors[i];
	std::sort(neighbors.begin(), neighbors.end());
	return neighbors;
}

cartesian_mesh make_mesh(index_t dims)
{
	return cartesian_mesh(dims, std::array<sindex_t, 3> { 0, 0, 0 }, std::array<sindex_t, 3> { 10, 10, 0 },
						  std::array<index_t, 3> { 2, 2, 1 });
}

} // namespace

TEST(UpdateCellNeighborsTest, NoAgentsDoesNotCrash)
{
	environment env(0.1, 2, 1, 1);
	auto mesh = make_mesh(2);
	kernels::openmp_solver::position_solver::update_cell_neighbors(env, mesh);
	EXPECT_EQ(env.get_agent_data().agents_count, 0);
}

TEST(UpdateCellNeighborsTest, SingleAgentHasNoNeighbors)
{
	environment env(0.1, 2, 1, 1);
	add_agent(env, { 0, 0 });

	auto mesh = make_mesh(2);
	kernels::openmp_solver::position_solver::update_cell_neighbors(env, mesh);

	ASSERT_EQ(env.get_agent_data().agents_count, 1);
	EXPECT_TRUE(env.get_agent_data().state_data.neighbors[0].empty());
}

TEST(UpdateCellNeighborsTest, DistanceEqualThresholdCountsAsNeighbor)
{
	environment env(0.1, 2, 1, 1);
	add_agent(env, { 0, 0 }, 1, 1, 1);
	add_agent(env, { 2, 0 }, 1, 1, 1); // adhesion_distance = 1*1 + 1*1 = 2

	auto mesh = make_mesh(2);
	kernels::openmp_solver::position_solver::update_cell_neighbors(env, mesh);

	EXPECT_EQ(sorted_neighbors(env, 0), (std::vector<index_t> { 1 }));
	EXPECT_EQ(sorted_neighbors(env, 1), (std::vector<index_t> { 0 }));
}

TEST(UpdateCellNeighborsTest, DistanceAboveThresholdIsNotNeighbor)
{
	environment env(0.1, 2, 1, 1);
	add_agent(env, { 0, 0 }, 1, 1, 1);
	add_agent(env, { 2.0001, 0 }, 1, 1, 1);

	auto mesh = make_mesh(2);
	kernels::openmp_solver::position_solver::update_cell_neighbors(env, mesh);

	EXPECT_TRUE(env.get_agent_data().state_data.neighbors[0].empty());
	EXPECT_TRUE(env.get_agent_data().state_data.neighbors[1].empty());
}

TEST(UpdateCellNeighborsTest, ClearsPreviousNeighborsAndRespectsMovableFlag)
{
	environment env(0.1, 2, 1, 1);
	add_agent(env, { 0, 0 }, 1, 1, 1); // movable
	add_agent(env, { 1, 0 }, 1, 0, 1); // immovable but within threshold
	add_agent(env, { 10, 0 }, 1, 1, 1);

	auto& data = env.get_agent_data();
	data.state_data.neighbors[0] = { 2, 12345 }; // garbage to prove clear
	data.state_data.neighbors[1] = { 0, 2 };	 // should be cleared, then skipped (immovable)

	auto mesh = make_mesh(2);
	kernels::openmp_solver::position_solver::update_cell_neighbors(env, mesh);

	EXPECT_EQ(sorted_neighbors(env, 0), (std::vector<index_t> { 1 }));
	EXPECT_TRUE(env.get_agent_data().state_data.neighbors[1].empty());
	EXPECT_TRUE(env.get_agent_data().state_data.neighbors[2].empty());

	// Integration: neighbor list should drive a non-zero force.
	std::fill(data.velocity.begin(), data.velocity.end(), static_cast<real_t>(0));
	kernels::openmp_solver::position_solver::update_cell_forces(env);
	EXPECT_NEAR(data.velocity[0], -data.velocity[2], 1e-6);
}
