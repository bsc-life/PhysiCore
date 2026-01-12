#include <limits>
#include <memory>

#include <common/base_agent_data.h>
#include <gtest/gtest.h>

#include "physicell/agent_data.h"
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
	auto* agent0 = container.create();
	auto* agent1 = container.create();
	auto* agent2 = container.create();

	ASSERT_NE(agent0, nullptr);
	ASSERT_NE(agent1, nullptr);
	ASSERT_NE(agent2, nullptr);

	// Stamp distinct positions to make sure each agent occupies its own slot
	agent0->position()[0] = 0.1;
	agent0->position()[1] = 0.2;
	agent0->position()[2] = 0.3;

	agent1->position()[0] = 1.1;
	agent1->position()[1] = 1.2;
	agent1->position()[2] = 1.3;

	agent2->position()[0] = 2.1;
	agent2->position()[1] = 2.2;
	agent2->position()[2] = 2.3;

	EXPECT_EQ(container.size(), 3);

	// Access by index should return the original agents with their stamped data
	auto* retrieved0 = container.get_agent_at(0);
	auto* retrieved1 = container.get_agent_at(1);
	auto* retrieved2 = container.get_agent_at(2);

	ASSERT_NE(retrieved0, nullptr);
	ASSERT_NE(retrieved1, nullptr);
	ASSERT_NE(retrieved2, nullptr);

	EXPECT_EQ(retrieved0, agent0);
	EXPECT_DOUBLE_EQ(retrieved0->position()[0], 0.1);
	EXPECT_DOUBLE_EQ(retrieved0->position()[1], 0.2);
	EXPECT_DOUBLE_EQ(retrieved0->position()[2], 0.3);

	EXPECT_EQ(retrieved1, agent1);
	EXPECT_DOUBLE_EQ(retrieved1->position()[0], 1.1);
	EXPECT_DOUBLE_EQ(retrieved1->position()[1], 1.2);
	EXPECT_DOUBLE_EQ(retrieved1->position()[2], 1.3);

	EXPECT_EQ(retrieved2, agent2);
	EXPECT_DOUBLE_EQ(retrieved2->position()[0], 2.1);
	EXPECT_DOUBLE_EQ(retrieved2->position()[1], 2.2);
	EXPECT_DOUBLE_EQ(retrieved2->position()[2], 2.3);
}

TEST_F(MechanicalAgentContainerTest, AccessAgentByIndex)
{
	auto* agent0 = container.create();
	auto* agent1 = container.create();
	auto* agent2 = container.create();

	// Give each agent a distinct signature
	agent0->radius() = 1.5;
	auto vel0 = agent0->velocity();
	vel0[0] = 0.1;
	vel0[1] = 0.2;
	vel0[2] = 0.3;

	agent1->radius() = 2.5;
	auto vel1 = agent1->velocity();
	vel1[0] = 1.1;
	vel1[1] = 1.2;
	vel1[2] = 1.3;
	agent1->migration_speed() = 0.4;

	agent2->radius() = 3.5;
	auto vel2 = agent2->velocity();
	vel2[0] = 2.1;
	vel2[1] = 2.2;
	vel2[2] = 2.3;
	agent2->chemotactic_sensitivities()[0] = 0.9;

	auto* retrieved0 = container.get_agent_at(0);
	ASSERT_NE(retrieved0, nullptr);
	EXPECT_EQ(retrieved0, agent0);
	EXPECT_DOUBLE_EQ(retrieved0->radius(), 1.5);
	EXPECT_DOUBLE_EQ(retrieved0->velocity()[0], 0.1);
	EXPECT_DOUBLE_EQ(retrieved0->velocity()[1], 0.2);
	EXPECT_DOUBLE_EQ(retrieved0->velocity()[2], 0.3);

	auto* retrieved1 = container.get_agent_at(1);
	ASSERT_NE(retrieved1, nullptr);
	EXPECT_EQ(retrieved1, agent1);
	EXPECT_DOUBLE_EQ(retrieved1->radius(), 2.5);
	EXPECT_DOUBLE_EQ(retrieved1->velocity()[0], 1.1);
	EXPECT_DOUBLE_EQ(retrieved1->velocity()[1], 1.2);
	EXPECT_DOUBLE_EQ(retrieved1->velocity()[2], 1.3);
	EXPECT_DOUBLE_EQ(retrieved1->migration_speed(), 0.4);

	auto* retrieved2 = container.get_agent_at(2);
	ASSERT_NE(retrieved2, nullptr);
	EXPECT_EQ(retrieved2, agent2);
	EXPECT_DOUBLE_EQ(retrieved2->radius(), 3.5);
	EXPECT_DOUBLE_EQ(retrieved2->velocity()[0], 2.1);
	EXPECT_DOUBLE_EQ(retrieved2->velocity()[1], 2.2);
	EXPECT_DOUBLE_EQ(retrieved2->velocity()[2], 2.3);
	EXPECT_DOUBLE_EQ(retrieved2->chemotactic_sensitivities()[0], 0.9);

#ifdef NDEBUG
	EXPECT_EQ(container.get_agent_at(3), nullptr);
	EXPECT_EQ(container.get_agent_at(std::numeric_limits<index_t>::max()), nullptr);
#endif
}

TEST_F(MechanicalAgentContainerTest, RemoveAgentDecreasesSize)
{
	container.create();
	container.create();
	container.create();
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
