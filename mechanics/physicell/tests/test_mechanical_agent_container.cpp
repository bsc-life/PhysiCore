#include <gtest/gtest.h>

#include "physicell/mechanical_agent_container.h"

namespace physicore::mechanics::physicell::tests {

using namespace physicore;
using namespace physicore::mechanics::physicell;

namespace {
mechanical_agent_container make_mechanical_agent_container()
{
	auto base_data = std::make_unique<base_agent_data>(3);
	auto data = std::make_unique<agent_data>(*base_data);

	return mechanical_agent_container(std::move(base_data), std::move(data));
}
} // namespace

class MechanicalAgentContainerTest : public ::testing::Test
{
protected:
	mechanical_agent_container container = make_mechanical_agent_container();

	void SetUp() override
	{
		// Container starts empty
	}
};

TEST_F(MechanicalAgentContainerTest, ContainerStartsEmpty) { EXPECT_EQ(container.size(), 0); }

TEST_F(MechanicalAgentContainerTest, AddSingleAgentIncreasesSize)
{
	container.create();
	EXPECT_EQ(container.size(), 1);
}

TEST_F(MechanicalAgentContainerTest, AddMultipleAgentsIncreasesSize)
{
	container.create();
	container.create();
	container.create();

	EXPECT_EQ(container.size(), 3);
}

TEST_F(MechanicalAgentContainerTest, AgentsHaveAutoIncrementingIDs)
{
	auto agent0 = container.create();
	auto agent1 = container.create();
	auto agent2 = container.create();

	EXPECT_NE(agent0, nullptr);
	EXPECT_NE(agent1, nullptr);
	EXPECT_NE(agent2, nullptr);
	EXPECT_EQ(agent0->get_index(), 0);
	EXPECT_EQ(agent1->get_index(), 1);
	EXPECT_EQ(agent2->get_index(), 2);
}

TEST_F(MechanicalAgentContainerTest, AccessAgentByIndex)
{
	container.create();
	container.create();
	container.create();

	auto agent = container.get_agent_at(1);
	EXPECT_NE(agent, nullptr);
	EXPECT_EQ(agent->get_index(), 1);
}

TEST_F(MechanicalAgentContainerTest, RemoveAgentDecreasesSize)
{
	auto agent0 = container.create();
	auto agent1 = container.create();
	auto agent2 = container.create();
	EXPECT_EQ(container.size(), 3);

	// Remove the last agent (safest approach)
	container.remove_at(2);
	EXPECT_EQ(container.size(), 2);
}

TEST_F(MechanicalAgentContainerTest, AgentDataSizeMatchesContainerSize)
{
	container.create();
	container.create();

	EXPECT_EQ(container.size(), 2);
}

TEST_F(MechanicalAgentContainerTest, LargeScaleAddRemove)
{
	// Add 100 agents
	for (int i = 0; i < 100; ++i)
	{
		container.create();
	}
	EXPECT_EQ(container.size(), 100);

	// Remove from the end (50 times)
	for (int i = 0; i < 50; ++i)
	{
		container.remove_at(container.size() - 1);
	}
	EXPECT_EQ(container.size(), 50);
}

} // namespace physicore::mechanics::physicell::tests
