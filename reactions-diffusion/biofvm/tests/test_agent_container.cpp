#include <gtest/gtest.h>

#include "agent.h"
#include "agent_container.h"

using namespace physicore;
using namespace physicore::biofvm;

agent_container make_agent_container()
{
	auto base_data = std::make_unique<base_agent_data>(3);
	auto data = std::make_unique<agent_data>(*base_data, 1);

	return agent_container(std::move(base_data), std::move(data));
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
	// After removal, agent1 pointer is invalid, but agent2 should still be valid
	EXPECT_NE(container.create(), nullptr);
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
	agent0->volume() = 1.0;
	agent0->secretion_rates()[0] = 0.1;
	agent0->saturation_densities()[0] = 0.2;
	agent0->uptake_rates()[0] = 0.3;
	agent0->net_export_rates()[0] = 0.4;
	agent0->internalized_substrates()[0] = 0.5;
	agent0->fraction_released_at_death()[0] = 0.6;
	agent0->fraction_transferred_when_ingested()[0] = 0.7;
	agent0->position()[0] = 0.8;

	agent1->volume() = 2.0;
	agent1->secretion_rates()[0] = 1.1;
	agent1->saturation_densities()[0] = 1.2;
	agent1->uptake_rates()[0] = 1.3;
	agent1->net_export_rates()[0] = 1.4;
	agent1->internalized_substrates()[0] = 1.5;
	agent1->fraction_released_at_death()[0] = 1.6;
	agent1->fraction_transferred_when_ingested()[0] = 1.7;
	agent1->position()[0] = 1.8;

	agent2->volume() = 3.0;
	agent2->secretion_rates()[0] = 2.1;
	agent2->saturation_densities()[0] = 2.2;
	agent2->uptake_rates()[0] = 2.3;
	agent2->net_export_rates()[0] = 2.4;
	agent2->internalized_substrates()[0] = 2.5;
	agent2->fraction_released_at_death()[0] = 2.6;
	agent2->fraction_transferred_when_ingested()[0] = 2.7;
	agent2->position()[0] = 2.8;

	int remove_idx = std::get<0>(GetParam());
	bool remove_via_pointer = std::get<1>(GetParam());

	if (remove_via_pointer)
		container.remove_agent(container.get_agent_at(remove_idx));
	else
		container.remove_at(remove_idx);

	// Verify remaining agents maintain their values
	if (remove_idx != 0)
	{
		EXPECT_DOUBLE_EQ(agent0->volume(), 1.0);
		EXPECT_DOUBLE_EQ(agent0->secretion_rates()[0], 0.1);
		EXPECT_DOUBLE_EQ(agent0->saturation_densities()[0], 0.2);
		EXPECT_DOUBLE_EQ(agent0->uptake_rates()[0], 0.3);
		EXPECT_DOUBLE_EQ(agent0->net_export_rates()[0], 0.4);
		EXPECT_DOUBLE_EQ(agent0->internalized_substrates()[0], 0.5);
		EXPECT_DOUBLE_EQ(agent0->fraction_released_at_death()[0], 0.6);
		EXPECT_DOUBLE_EQ(agent0->fraction_transferred_when_ingested()[0], 0.7);
		EXPECT_DOUBLE_EQ(agent0->position()[0], 0.8);
	}
	if (remove_idx != 1)
	{
		EXPECT_DOUBLE_EQ(agent1->volume(), 2.0);
		EXPECT_DOUBLE_EQ(agent1->secretion_rates()[0], 1.1);
		EXPECT_DOUBLE_EQ(agent1->saturation_densities()[0], 1.2);
		EXPECT_DOUBLE_EQ(agent1->uptake_rates()[0], 1.3);
		EXPECT_DOUBLE_EQ(agent1->net_export_rates()[0], 1.4);
		EXPECT_DOUBLE_EQ(agent1->internalized_substrates()[0], 1.5);
		EXPECT_DOUBLE_EQ(agent1->fraction_released_at_death()[0], 1.6);
		EXPECT_DOUBLE_EQ(agent1->fraction_transferred_when_ingested()[0], 1.7);
		EXPECT_DOUBLE_EQ(agent1->position()[0], 1.8);
	}
	if (remove_idx != 2)
	{
		EXPECT_DOUBLE_EQ(agent2->volume(), 3.0);
		EXPECT_DOUBLE_EQ(agent2->secretion_rates()[0], 2.1);
		EXPECT_DOUBLE_EQ(agent2->saturation_densities()[0], 2.2);
		EXPECT_DOUBLE_EQ(agent2->uptake_rates()[0], 2.3);
		EXPECT_DOUBLE_EQ(agent2->net_export_rates()[0], 2.4);
		EXPECT_DOUBLE_EQ(agent2->internalized_substrates()[0], 2.5);
		EXPECT_DOUBLE_EQ(agent2->fraction_released_at_death()[0], 2.6);
		EXPECT_DOUBLE_EQ(agent2->fraction_transferred_when_ingested()[0], 2.7);
		EXPECT_DOUBLE_EQ(agent2->position()[0], 2.8);
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
	agent0->volume() = 1.5;
	agent0->secretion_rates()[0] = 0.1;

	agent1->volume() = 2.5;
	agent1->secretion_rates()[0] = 0.2;

	agent2->volume() = 3.5;
	agent2->secretion_rates()[0] = 0.3;

	// Test get_agent_at for each position
	auto* retrieved0 = container.get_agent_at(0);
	ASSERT_NE(retrieved0, nullptr);
	EXPECT_DOUBLE_EQ(retrieved0->volume(), 1.5);
	EXPECT_DOUBLE_EQ(retrieved0->secretion_rates()[0], 0.1);
	EXPECT_EQ(retrieved0, agent0);

	auto* retrieved1 = container.get_agent_at(1);
	ASSERT_NE(retrieved1, nullptr);
	EXPECT_DOUBLE_EQ(retrieved1->volume(), 2.5);
	EXPECT_DOUBLE_EQ(retrieved1->secretion_rates()[0], 0.2);
	EXPECT_EQ(retrieved1, agent1);

	auto* retrieved2 = container.get_agent_at(2);
	ASSERT_NE(retrieved2, nullptr);
	EXPECT_DOUBLE_EQ(retrieved2->volume(), 3.5);
	EXPECT_DOUBLE_EQ(retrieved2->secretion_rates()[0], 0.3);
	EXPECT_EQ(retrieved2, agent2);

#ifdef NDEBUG
	// Test out of bounds access
	EXPECT_EQ(container.get_agent_at(3), nullptr);
	EXPECT_EQ(container.get_agent_at(std::numeric_limits<index_t>::max()), nullptr);
#endif
}
