#include <gtest/gtest.h>

#include "base_agent_container.h"

using namespace physicore;

TEST(BaseAgentContainerTest, CreateAndRemove)
{
	base_agent_container container(std::make_unique<base_agent_data>());
	auto* agent1 = container.create();
	EXPECT_EQ(container.size(), 1);
	EXPECT_NE(agent1, nullptr);
	auto* agent2 = container.create();
	EXPECT_EQ(container.size(), 2);
	EXPECT_NE(agent2, nullptr);
	EXPECT_NE(agent1, agent2);
	container.remove_agent(agent1);
	EXPECT_EQ(container.size(), 1);
	// After removal, agent1 pointer is invalid, but agent2 should still be valid
	EXPECT_NE(container.create(), nullptr);
}

class RemoveAgentTest : public ::testing::TestWithParam<int>
{};

TEST_P(RemoveAgentTest, RemoveAgentsAndCheckPositions)
{
	base_agent_container container(std::make_unique<base_agent_data>());
	auto* agent0 = container.create();
	auto* agent1 = container.create();
	auto* agent2 = container.create();

	agent0->position()[0] = 1.0;
	agent0->position()[1] = 2.0; // agent 0
	agent1->position()[0] = 3.0;
	agent1->position()[1] = 4.0; // agent 1
	agent2->position()[0] = 5.0;
	agent2->position()[1] = 6.0; // agent 2

	const int remove_idx = GetParam();
	container.remove_at(remove_idx);

	if (remove_idx != 0)
	{
		EXPECT_EQ(agent0->position()[0], 1.0);
		EXPECT_EQ(agent0->position()[1], 2.0);
	}
	if (remove_idx != 1)
	{
		EXPECT_EQ(agent1->position()[0], 3.0);
		EXPECT_EQ(agent1->position()[1], 4.0);
	}
	if (remove_idx != 2)
	{
		EXPECT_EQ(agent2->position()[0], 5.0);
		EXPECT_EQ(agent2->position()[1], 6.0);
	}
}

INSTANTIATE_TEST_SUITE_P(BaseAgentContainerTest, RemoveAgentTest, ::testing::Values(0, 1, 2));
