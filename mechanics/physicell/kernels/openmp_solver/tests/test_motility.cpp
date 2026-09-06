#include <gtest/gtest.h>
#include <physicell/openmp_solver/position_solver.h>

using namespace physicore;
using namespace physicore::mechanics::physicell;

class agent_retriever : public generic_agent_solver<mechanical_agent>
{
public:
	using generic_agent_solver<mechanical_agent>::retrieve_agent_data;
};

class UpdateMotilityTest : public ::testing::TestWithParam<index_t>
{};

TEST_P(UpdateMotilityTest, AddsMotilityVectorToVelocity)
{
	const index_t dims = GetParam();
	environment env({ dims, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 } }, 1, 1, 0);
	auto* agent = env.agents->create();
	agent->is_motile() = 1;

	for (index_t d = 0; d < dims; ++d)
	{
		agent->motility_vector()[d] = static_cast<real_t>(d + 1);
		agent->velocity()[d] = static_cast<real_t>(10 + d);
	}

	kernels::openmp_solver::position_solver solver;
	solver.update_motility(env);

	for (index_t d = 0; d < dims; ++d)
		EXPECT_FLOAT_EQ(agent->velocity()[d], static_cast<real_t>(11 + 2 * d));
}

TEST_P(UpdateMotilityTest, NonMotileAgentIsUnchanged)
{
	const index_t dims = GetParam();
	environment env({ dims, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 } }, 1, 1, 0);
	auto* agent = env.agents->create();

	for (index_t d = 0; d < dims; ++d)
	{
		agent->motility_vector()[d] = static_cast<real_t>(d + 1);
		agent->velocity()[d] = static_cast<real_t>(10 + d);
	}

	kernels::openmp_solver::position_solver solver;
	solver.update_motility(env);

	for (index_t d = 0; d < dims; ++d)
	{
		EXPECT_FLOAT_EQ(agent->motility_vector()[d], static_cast<real_t>(d + 1));
		EXPECT_FLOAT_EQ(agent->velocity()[d], static_cast<real_t>(10 + d));
	}
}

TEST(UpdateMotilityTest, RefreshesBiasedDirectionAndCallsCallback)
{
	const index_t dims = 2;
	environment env({ dims, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 } }, 1, 1, 2);
	auto* agent = env.agents->create();
	agent->is_motile() = 1;
	agent->persistence_time() = 1;
	agent->migration_speed() = 10;
	agent->migration_bias() = 1;
	agent->agent_type_index() = 0;

	bool callback_called = false;
	index_t callback_type = -1;
	auto& data = agent_retriever().retrieve_agent_data(*env.agents);
	data.motility_data.direction_update_funcs[0] = [&](index_t type) {
		callback_called = true;
		callback_type = type;
		agent->migration_bias_direction()[0] = 3;
		agent->migration_bias_direction()[1] = 4;
	};

	kernels::openmp_solver::position_solver solver;
	solver.update_motility(env);

	EXPECT_TRUE(callback_called);
	EXPECT_EQ(callback_type, 0);
	EXPECT_FLOAT_EQ(agent->motility_vector()[0], 6);
	EXPECT_FLOAT_EQ(agent->motility_vector()[1], 8);
	EXPECT_FLOAT_EQ(agent->velocity()[0], 6);
	EXPECT_FLOAT_EQ(agent->velocity()[1], 8);
}

INSTANTIATE_TEST_SUITE_P(Dimensions, UpdateMotilityTest, ::testing::Values(index_t(1), index_t(2), index_t(3)));
