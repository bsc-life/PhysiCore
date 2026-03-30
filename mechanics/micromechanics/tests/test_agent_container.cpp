#include <memory>

#include <gtest/gtest.h>

#include "micromechanics/agent_container.h"
#include "micromechanics/agent_data.h"

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

static agent_container make_agent_container()
{
	auto base_data = std::make_unique<base_agent_data>(3);
	auto mech_data = std::make_unique<agent_data>(*base_data);
	return agent_container(std::move(base_data), std::move(mech_data));
}

TEST(AgentContainerTest, CreateIncreasesSize)
{
	agent_container container = make_agent_container();

	auto* a0 = container.create();
	ASSERT_NE(a0, nullptr);
	EXPECT_EQ(container.size(), 1U);

	auto* a1 = container.create();
	ASSERT_NE(a1, nullptr);
	EXPECT_EQ(container.size(), 2U);
}

TEST(AgentContainerTest, CreateAndRemove)
{
	agent_container container = make_agent_container();
	agent* agent1 = container.create();
	EXPECT_NE(agent1, nullptr);
	agent* agent2 = container.create();
	EXPECT_NE(agent2, nullptr);
	EXPECT_NE(agent1, agent2);
	container.remove_agent(agent1);
	EXPECT_EQ(container.size(), 1U);
	// After removal, agent1 pointer is invalid, but we can create a new agent
	EXPECT_NE(container.create(), nullptr);
	EXPECT_EQ(container.size(), 2U);
}

class RemoveAgentTest : public ::testing::TestWithParam<std::tuple<int, bool>>
{};

TEST_P(RemoveAgentTest, RemoveAgentsAndCheckProperties)
{
	agent_container container = make_agent_container();

	// Create test agents
	auto* agent0 = container.create();
	auto* agent1 = container.create();
	auto* agent2 = container.create();

	// Set unique values for each agent to verify correct data management
	agent0->compartment_type() = 1;
	agent0->cell_id() = 10;
	agent0->position()[0] = 0.8;
	agent0->velocity()[0] = 0.9;

	agent1->compartment_type() = 2;
	agent1->cell_id() = 20;
	agent1->position()[0] = 1.8;
	agent1->velocity()[0] = 1.9;

	agent2->compartment_type() = 3;
	agent2->cell_id() = 30;
	agent2->position()[0] = 2.8;
	agent2->velocity()[0] = 2.9;

	int const remove_idx = std::get<0>(GetParam());
	bool const remove_via_pointer = std::get<1>(GetParam());

	if (remove_via_pointer)
		container.remove_agent(container.get_agent_at(remove_idx));
	else
		container.remove_at(remove_idx);

	// Verify remaining agents maintain their values
	if (remove_idx != 0)
	{
		EXPECT_EQ(agent0->compartment_type(), 1);
		EXPECT_EQ(agent0->cell_id(), 10);
		EXPECT_DOUBLE_EQ(agent0->position()[0], 0.8);
		EXPECT_DOUBLE_EQ(agent0->velocity()[0], 0.9);
	}
	if (remove_idx != 1)
	{
		EXPECT_EQ(agent1->compartment_type(), 2);
		EXPECT_EQ(agent1->cell_id(), 20);
		EXPECT_DOUBLE_EQ(agent1->position()[0], 1.8);
		EXPECT_DOUBLE_EQ(agent1->velocity()[0], 1.9);
	}
	if (remove_idx != 2)
	{
		EXPECT_EQ(agent2->compartment_type(), 3);
		EXPECT_EQ(agent2->cell_id(), 30);
		EXPECT_DOUBLE_EQ(agent2->position()[0], 2.8);
		EXPECT_DOUBLE_EQ(agent2->velocity()[0], 2.9);
	}
}

INSTANTIATE_TEST_SUITE_P(AgentContainerTest, RemoveAgentTest,
						 ::testing::Combine(::testing::Values(0, 1, 2), ::testing::Values(true, false)));

TEST(AgentContainerTest, GetAgentAt)
{
	agent_container container = make_agent_container();

	// Create multiple agents and set unique values
	auto* agent0 = container.create();
	auto* agent1 = container.create();
	auto* agent2 = container.create();

	// Set distinctive values for each agent
	agent0->compartment_type() = 11;
	agent0->cell_id() = 101;

	agent1->compartment_type() = 22;
	agent1->cell_id() = 202;

	agent2->compartment_type() = 33;
	agent2->cell_id() = 303;

	// Test get_agent_at for each position
	auto* retrieved0 = container.get_agent_at(0);
	ASSERT_NE(retrieved0, nullptr);
	EXPECT_EQ(retrieved0->compartment_type(), 11);
	EXPECT_EQ(retrieved0->cell_id(), 101);
	EXPECT_EQ(retrieved0, agent0);

	auto* retrieved1 = container.get_agent_at(1);
	ASSERT_NE(retrieved1, nullptr);
	EXPECT_EQ(retrieved1->compartment_type(), 22);
	EXPECT_EQ(retrieved1->cell_id(), 202);
	EXPECT_EQ(retrieved1, agent1);

	auto* retrieved2 = container.get_agent_at(2);
	ASSERT_NE(retrieved2, nullptr);
	EXPECT_EQ(retrieved2->compartment_type(), 33);
	EXPECT_EQ(retrieved2->cell_id(), 303);
	EXPECT_EQ(retrieved2, agent2);

#ifdef NDEBUG
	// Test out of bounds access
	EXPECT_EQ(container.get_agent_at(3), nullptr);
	EXPECT_EQ(container.get_agent_at(-1), nullptr);
#endif
}
