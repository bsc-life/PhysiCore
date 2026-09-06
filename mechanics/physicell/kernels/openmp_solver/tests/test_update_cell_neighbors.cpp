#include <algorithm>

#include <gtest/gtest.h>
#include <physicell/openmp_solver/position_solver.h>

using namespace physicore::mechanics::physicell;
using physicore::cartesian_mesh;
using physicore::index_t;
using physicore::real_t;
using physicore::sindex_t;

namespace {

class agent_retriever : public physicore::generic_agent_solver<mechanical_agent>
{
public:
	using physicore::generic_agent_solver<mechanical_agent>::retrieve_agent_data;
};

mechanical_agent_data& retrieve_environment_agent_data(environment& env)
{
	if (env.agents == nullptr)
	{
		throw std::runtime_error("environment has no agents");
	}

	return agent_retriever().retrieve_agent_data(*env.agents);
}

kernels::openmp_solver::position_solver& position_solver_instance()
{
	static kernels::openmp_solver::position_solver solver;
	return solver;
}

void add_agent(environment& env, std::initializer_list<real_t> pos, real_t radius = 1, std::uint8_t movable = 1,
			   real_t rel_max_adhesion_dist = 1)
{
	ASSERT_TRUE(env.agents != nullptr);
	env.agents->create();

	auto& data = retrieve_environment_agent_data(env);
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
	auto neighbors = retrieve_environment_agent_data(env).state_data.neighbors[i];
	std::ranges::sort(neighbors);
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
	environment env(make_mesh(2), 1, 1, 0.1);
	auto mesh = make_mesh(2);
	position_solver_instance().update_cell_neighbors(env, mesh);
	EXPECT_EQ(retrieve_environment_agent_data(env).agents_count, 0);
}

TEST(UpdateCellNeighborsTest, SingleAgentHasNoNeighbors)
{
	environment env(make_mesh(2), 1, 1, 0.1);
	add_agent(env, { 0, 0 });

	auto mesh = make_mesh(2);
	position_solver_instance().update_cell_neighbors(env, mesh);

	ASSERT_EQ(retrieve_environment_agent_data(env).agents_count, 1);
	EXPECT_TRUE(retrieve_environment_agent_data(env).state_data.neighbors[0].empty());
}

TEST(UpdateCellNeighborsTest, DistanceEqualThresholdCountsAsNeighbor)
{
	environment env(make_mesh(2), 1, 1, 0.1);
	add_agent(env, { 0, 0 }, 1, 1, 1);
	add_agent(env, { 2, 0 }, 1, 1, 1); // adhesion_distance = 1*1 + 1*1 = 2

	auto mesh = make_mesh(2);
	position_solver_instance().update_cell_neighbors(env, mesh);

	EXPECT_EQ(sorted_neighbors(env, 0), (std::vector<index_t> { 1 }));
	EXPECT_EQ(sorted_neighbors(env, 1), (std::vector<index_t> { 0 }));
}

TEST(UpdateCellNeighborsTest, DistanceAboveThresholdIsNotNeighbor)
{
	environment env(make_mesh(2), 1, 1, 0.1);
	add_agent(env, { 0, 0 }, 1, 1, 1);
	add_agent(env, { 2.0001, 0 }, 1, 1, 1);

	auto mesh = make_mesh(2);
	position_solver_instance().update_cell_neighbors(env, mesh);

	EXPECT_TRUE(retrieve_environment_agent_data(env).state_data.neighbors[0].empty());
	EXPECT_TRUE(retrieve_environment_agent_data(env).state_data.neighbors[1].empty());
}

TEST(UpdateCellNeighborsTest, ClearsPreviousNeighborsAndRespectsMovableFlag)
{
	environment env(make_mesh(2), 1, 1, 0.1);
	add_agent(env, { 0, 0 }, 1, 1, 1); // movable
	add_agent(env, { 1, 0 }, 1, 0, 1); // immovable but within threshold
	add_agent(env, { 10, 0 }, 1, 1, 1);

	auto& data = retrieve_environment_agent_data(env);
	data.state_data.neighbors[0] = { 2, 12345 }; // garbage to prove clear
	data.state_data.neighbors[1] = { 0, 2 };	 // should be cleared, then skipped (immovable)

	auto mesh = make_mesh(2);
	position_solver_instance().update_cell_neighbors(env, mesh);

	EXPECT_EQ(sorted_neighbors(env, 0), (std::vector<index_t> { 1 }));
	EXPECT_TRUE(retrieve_environment_agent_data(env).state_data.neighbors[1].empty());
	EXPECT_TRUE(retrieve_environment_agent_data(env).state_data.neighbors[2].empty());

	// Integration: neighbor list should drive a non-zero force.
	std::ranges::fill(data.velocity, static_cast<real_t>(0));
	position_solver_instance().update_cell_forces(env);
	EXPECT_NEAR(data.velocity[0], -data.velocity[2], 1e-6);
}
