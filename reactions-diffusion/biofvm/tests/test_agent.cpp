#include <gtest/gtest.h>

#include "agent.h"
#include "agent_data.h"
#include "base_agent_data.h"

using namespace physicore;
using namespace physicore::biofvm;

class AgentTest : public ::testing::Test
{
protected:
	void SetUp() override
	{
		data.base_data.dims = 3;
		data.base_data.agents_count = 0;
		data.substrate_count = 2; // Testing with 2 substrates
		data.add_component<secretion_uptake_component>();
		data.add_component<net_export_component>();
		data.add(); // Add first agent
	}

	base_agent_data base_data;
	agent_data data = agent_data(base_data, 2);
};

TEST_F(AgentTest, SecretionRates)
{
	secretion_uptake_agent_mixin<agent> test_agent(0, data);
	auto rates = test_agent.secretion_rates();
	ASSERT_EQ(rates.size(), 2);

	// Test setting and getting values
	rates[0] = 1.5;
	rates[1] = 2.5;
	EXPECT_EQ(test_agent.secretion_rates()[0], 1.5);
	EXPECT_EQ(test_agent.secretion_rates()[1], 2.5);
}

TEST_F(AgentTest, SaturationDensities)
{
	secretion_uptake_agent_mixin<agent> test_agent(0, data);
	auto densities = test_agent.saturation_densities();
	ASSERT_EQ(densities.size(), 2);

	densities[0] = 10.0;
	densities[1] = 20.0;
	EXPECT_EQ(test_agent.saturation_densities()[0], 10.0);
	EXPECT_EQ(test_agent.saturation_densities()[1], 20.0);
}

TEST_F(AgentTest, UptakeRates)
{
	secretion_uptake_agent_mixin<agent> test_agent(0, data);
	auto rates = test_agent.uptake_rates();
	ASSERT_EQ(rates.size(), 2);

	rates[0] = 0.5;
	rates[1] = 1.0;
	EXPECT_EQ(test_agent.uptake_rates()[0], 0.5);
	EXPECT_EQ(test_agent.uptake_rates()[1], 1.0);
}

TEST_F(AgentTest, NetExportRates)
{
	net_export_agent_mixin<agent> test_agent(0, data);
	auto rates = test_agent.net_export_rates();
	ASSERT_EQ(rates.size(), 2);

	rates[0] = 3.0;
	rates[1] = 4.0;
	EXPECT_EQ(test_agent.net_export_rates()[0], 3.0);
	EXPECT_EQ(test_agent.net_export_rates()[1], 4.0);
}

TEST_F(AgentTest, InternalizedSubstrates)
{
	agent test_agent(0, data);
	auto substrates = test_agent.internalized_substrates();
	ASSERT_EQ(substrates.size(), 2);

	substrates[0] = 5.5;
	substrates[1] = 6.5;
	EXPECT_EQ(test_agent.internalized_substrates()[0], 5.5);
	EXPECT_EQ(test_agent.internalized_substrates()[1], 6.5);
}

TEST_F(AgentTest, FractionReleasedAtDeath)
{
	agent test_agent(0, data);
	auto fractions = test_agent.fraction_released_at_death();
	ASSERT_EQ(fractions.size(), 2);

	fractions[0] = 0.75;
	fractions[1] = 0.85;
	EXPECT_EQ(test_agent.fraction_released_at_death()[0], 0.75);
	EXPECT_EQ(test_agent.fraction_released_at_death()[1], 0.85);
}

TEST_F(AgentTest, FractionTransferredWhenIngested)
{
	agent test_agent(0, data);
	auto fractions = test_agent.fraction_transferred_when_ingested();
	ASSERT_EQ(fractions.size(), 2);

	fractions[0] = 0.25;
	fractions[1] = 0.35;
	EXPECT_EQ(test_agent.fraction_transferred_when_ingested()[0], 0.25);
	EXPECT_EQ(test_agent.fraction_transferred_when_ingested()[1], 0.35);
}

TEST_F(AgentTest, Volume)
{
	agent test_agent(0, data);
	test_agent.volume() = 100.5;
	EXPECT_EQ(test_agent.volume(), 100.5);
}

TEST_F(AgentTest, MultipleAgents)
{
	data.add(); // Add second agent
	secretion_uptake_agent_mixin<agent> agent0(0, data);
	secretion_uptake_agent_mixin<agent> agent1(1, data);

	// Test that the agents access different memory locations
	agent0.secretion_rates()[0] = 1.0;
	agent0.secretion_rates()[1] = 2.0;
	agent1.secretion_rates()[0] = 3.0;
	agent1.secretion_rates()[1] = 4.0;

	EXPECT_EQ(agent0.secretion_rates()[0], 1.0);
	EXPECT_EQ(agent0.secretion_rates()[1], 2.0);
	EXPECT_EQ(agent1.secretion_rates()[0], 3.0);
	EXPECT_EQ(agent1.secretion_rates()[1], 4.0);

	// Test volume for multiple agents
	agent0.volume() = 100.0;
	agent1.volume() = 200.0;
	EXPECT_EQ(agent0.volume(), 100.0);
	EXPECT_EQ(agent1.volume(), 200.0);
}
