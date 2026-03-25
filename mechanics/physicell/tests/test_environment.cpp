#include <memory>
#include <stdexcept>

#include <common/generic_agent_solver.h>

#include <common/types.h>
#include <gtest/gtest.h>

#include "environment.h"

using namespace physicore;
using namespace physicore::mechanics::physicell;

namespace {

class agent_retriever : public generic_agent_solver<mechanical_agent>
{};

mechanical_agent_data& retrieve_environment_agent_data(environment& env)
{
	if (env.agents == nullptr)
	{
		throw std::runtime_error("environment has no agents");
	}

	return agent_retriever().retrieve_agent_data(*env.agents);
}

constexpr real_t test_timestep = 0.1;
constexpr index_t test_dims = 2;
constexpr index_t test_agent_types = 1;
constexpr index_t test_substrates = 1;

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

TEST(EnvironmentTest, RunSingleTimestepUsesSolverWhenProvided)
{
	environment env(test_timestep, test_dims, test_agent_types, test_substrates);
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
	environment env(test_timestep, test_dims, test_agent_types, test_substrates);
	EXPECT_NO_THROW(env.serialize_state(2.5));
}

TEST(EnvironmentTest, SerializeStateUsesSerializerWhenProvided)
{
	environment env(test_timestep, test_dims, test_agent_types, test_substrates);
	auto serializer = std::make_unique<TestSerializer>();
	auto* serializer_ptr = serializer.get();
	env.serializer = std::move(serializer);

	env.serialize_state(3.25);

	EXPECT_EQ(serializer_ptr->calls, 1);
	EXPECT_DOUBLE_EQ(serializer_ptr->last_time, 3.25);
}

TEST(EnvironmentTest, RetrieveAgentDataReturnsContainerData)
{
	environment env(0.1, 2, 4, 3);
	auto& data = retrieve_environment_agent_data(env);

	ASSERT_NE(env.agents, nullptr);
	auto* expected = std::get<std::unique_ptr<mechanical_agent_data>>(env.agents->agent_datas).get();
	ASSERT_NE(expected, nullptr);
	EXPECT_EQ(&data, expected);
	EXPECT_EQ(data.agent_types_count, 4);
	EXPECT_EQ(data.substrates_count, 3);
	EXPECT_EQ(data.base_data.dims, 2);
}

TEST(EnvironmentTest, RetrieveAgentDataThrowsWhenAgentsMissing)
{
	environment env(test_timestep, test_dims, test_agent_types, test_substrates);
	env.agents.reset();

	EXPECT_THROW(retrieve_environment_agent_data(env), std::runtime_error);
}
