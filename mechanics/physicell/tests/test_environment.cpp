#include <memory>

#include <gtest/gtest.h>

#include <common/types.h>

#include "environment.h"

using namespace physicore;
using namespace physicore::mechanics::physicell;

namespace {

class TestSolver final : public solver
{
public:
	void initialize(environment& /*e*/) override {}

	void solve(environment& e, index_t iterations) override
	{
		last_env = &e;
		last_iterations = iterations;
		++calls;
	}

	int calls = 0;
	environment* last_env = nullptr;
	index_t last_iterations = 0;
};

class TestSerializer final : public serializer
{
public:
	void serialize(real_t current_time) override
	{
		last_time = current_time;
		++calls;
	}

	int calls = 0;
	real_t last_time = 0.0;
};

} // namespace

TEST(EnvironmentTest, RunSingleTimestepWithoutSolver)
{
	environment env(0.1);
	EXPECT_NO_THROW(env.run_single_timestep());
}

TEST(EnvironmentTest, RunSingleTimestepUsesSolverWhenProvided)
{
	environment env(0.1);
	auto solver = std::make_unique<TestSolver>();
	auto* solver_ptr = solver.get();
	env.solver = std::move(solver);

	env.run_single_timestep();

	EXPECT_EQ(solver_ptr->calls, 1);
	EXPECT_EQ(solver_ptr->last_env, &env);
	EXPECT_EQ(solver_ptr->last_iterations, 1);
}

TEST(EnvironmentTest, SerializeStateWithoutSerializer)
{
	environment env(0.1);
	EXPECT_NO_THROW(env.serialize_state(2.5));
}

TEST(EnvironmentTest, SerializeStateUsesSerializerWhenProvided)
{
	environment env(0.1);
	auto serializer = std::make_unique<TestSerializer>();
	auto* serializer_ptr = serializer.get();
	env.serializer = std::move(serializer);

	env.serialize_state(3.25);

	EXPECT_EQ(serializer_ptr->calls, 1);
	EXPECT_DOUBLE_EQ(serializer_ptr->last_time, 3.25);
}

TEST(EnvironmentTest, GetAgentDataReturnsContainerData)
{
	environment env(0.1, 2, 4, 3);
	auto& data = env.get_agent_data();

	ASSERT_NE(env.agents, nullptr);
	auto* expected = std::get<std::unique_ptr<mechanical_agent_data>>(env.agents->agent_datas).get();
	ASSERT_NE(expected, nullptr);
	EXPECT_EQ(&data, expected);
	EXPECT_EQ(data.agent_types_count, 4);
	EXPECT_EQ(data.substrates_count, 3);
	EXPECT_EQ(data.base_data.dims, 2);
}

TEST(EnvironmentTest, GetAgentDataThrowsWhenAgentsMissing)
{
	environment env(0.1);
	env.agents.reset();

	EXPECT_THROW(env.get_agent_data(), std::runtime_error);
}
