#include <common/base_agent_data.h>
#include <gtest/gtest.h>

#include "agent_data.h"

using namespace physicore::biofvm;
using namespace physicore;

// Helper to create a base_agent_data with N agents
base_agent_data make_base_agent_data(index_t count)
{
	base_agent_data data;
	for (index_t i = 0; i < count; ++i)
	{
		data.add();
	}
	return data;
}

TEST(AgentDataTest, AddInitializesVectorsCorrectly)
{
	base_agent_data base = make_base_agent_data(0);
	index_t substrate_count = 3;
	agent_data data(base, substrate_count);

	data.add(); // adds one agent

	EXPECT_EQ(data.agents_count, 1);
	EXPECT_EQ(data.secretion_rates.size(), substrate_count);
	EXPECT_EQ(data.saturation_densities.size(), substrate_count);
	EXPECT_EQ(data.uptake_rates.size(), substrate_count);
	EXPECT_EQ(data.net_export_rates.size(), substrate_count);
	EXPECT_EQ(data.internalized_substrates.size(), substrate_count);
	EXPECT_EQ(data.fraction_released_at_death.size(), substrate_count);
	EXPECT_EQ(data.fraction_transferred_when_ingested.size(), substrate_count);
	EXPECT_EQ(data.volumes.size(), 1);

	data.add(); // adds another agent

	EXPECT_EQ(data.agents_count, 2);
	EXPECT_EQ(data.secretion_rates.size(), substrate_count * 2);
	EXPECT_EQ(data.saturation_densities.size(), substrate_count * 2);
	EXPECT_EQ(data.uptake_rates.size(), substrate_count * 2);
	EXPECT_EQ(data.net_export_rates.size(), substrate_count * 2);
	EXPECT_EQ(data.internalized_substrates.size(), substrate_count * 2);
	EXPECT_EQ(data.fraction_released_at_death.size(), substrate_count * 2);
	EXPECT_EQ(data.fraction_transferred_when_ingested.size(), substrate_count * 2);
	EXPECT_EQ(data.volumes.size(), 2);
}

TEST(AgentDataTest, RemoveShrinksVectorsCorrectly)
{
	base_agent_data base = make_base_agent_data(0);
	index_t substrate_count = 2;
	agent_data data(base, substrate_count);

	data.add(); // agent 0
	data.add(); // agent 1
	data.add(); // agent 2

	index_t expected_size = 3;

	EXPECT_EQ(data.agents_count, expected_size);
	EXPECT_EQ(data.secretion_rates.size(), substrate_count * expected_size);
	EXPECT_EQ(data.saturation_densities.size(), substrate_count * expected_size);
	EXPECT_EQ(data.uptake_rates.size(), substrate_count * expected_size);
	EXPECT_EQ(data.net_export_rates.size(), substrate_count * expected_size);
	EXPECT_EQ(data.internalized_substrates.size(), substrate_count * expected_size);
	EXPECT_EQ(data.fraction_released_at_death.size(), substrate_count * expected_size);
	EXPECT_EQ(data.fraction_transferred_when_ingested.size(), substrate_count * expected_size);
	EXPECT_EQ(data.volumes.size(), expected_size);

	// Remove agent at position 1
	data.remove_at(1);

	expected_size = 2;

	EXPECT_EQ(data.agents_count, expected_size);
	EXPECT_EQ(data.secretion_rates.size(), substrate_count * expected_size);
	EXPECT_EQ(data.saturation_densities.size(), substrate_count * expected_size);
	EXPECT_EQ(data.uptake_rates.size(), substrate_count * expected_size);
	EXPECT_EQ(data.net_export_rates.size(), substrate_count * expected_size);
	EXPECT_EQ(data.internalized_substrates.size(), substrate_count * expected_size);
	EXPECT_EQ(data.fraction_released_at_death.size(), substrate_count * expected_size);
	EXPECT_EQ(data.fraction_transferred_when_ingested.size(), substrate_count * expected_size);
	EXPECT_EQ(data.volumes.size(), expected_size);

	// Remove agent at position 0
	data.remove_at(0);

	expected_size = 1;

	EXPECT_EQ(data.agents_count, expected_size);
	EXPECT_EQ(data.secretion_rates.size(), substrate_count * expected_size);
	EXPECT_EQ(data.saturation_densities.size(), substrate_count * expected_size);
	EXPECT_EQ(data.uptake_rates.size(), substrate_count * expected_size);
	EXPECT_EQ(data.net_export_rates.size(), substrate_count * expected_size);
	EXPECT_EQ(data.internalized_substrates.size(), substrate_count * expected_size);
	EXPECT_EQ(data.fraction_released_at_death.size(), substrate_count * expected_size);
	EXPECT_EQ(data.fraction_transferred_when_ingested.size(), substrate_count * expected_size);
	EXPECT_EQ(data.volumes.size(), expected_size);
}
