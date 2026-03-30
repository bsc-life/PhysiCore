#include <common/base_agent_data.h>
#include <gtest/gtest.h>

#include "micromechanics/agent_data.h"

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

class AgentDataTest : public ::testing::Test
{
protected:
	void SetUp() override
	{
		base_data.dims = 3;
		base_data.agents_count = 0;
	}

	base_agent_data base_data;
};

TEST_F(AgentDataTest, AddIncreasesCount)
{
	agent_data data(base_data);
	EXPECT_EQ(data.agents_count, 0);

	data.add();
	EXPECT_EQ(data.agents_count, 1);

	data.add();
	EXPECT_EQ(data.agents_count, 2);

	data.add();
	EXPECT_EQ(data.agents_count, 3);
}

TEST_F(AgentDataTest, AddResizesAllVectors)
{
	agent_data data(base_data);
	data.add();

	// Scalar fields
	EXPECT_EQ(data.compartment_types.size(), 1U);
	EXPECT_EQ(data.cell_ids.size(), 1U);

	// Vector fields (3 components each)
	EXPECT_EQ(data.velocities.size(), 3U);
	EXPECT_EQ(data.previous_velocities.size(), 3U);
	EXPECT_EQ(data.forces.size(), 3U);

	// Nested vectors
	EXPECT_EQ(data.spring_attachments.size(), 1U);
}

TEST_F(AgentDataTest, AddInitializesDefaultValues)
{
	agent_data data(base_data);
	data.add();

	// Check default values
	EXPECT_EQ(data.compartment_types[0], 0);
	EXPECT_EQ(data.cell_ids[0], static_cast<index_t>(-1)); // Standalone agent
	EXPECT_DOUBLE_EQ(data.velocities[0], 0.0);
	EXPECT_DOUBLE_EQ(data.forces[0], 0.0);
}

TEST_F(AgentDataTest, RemoveAtDecreasesCount)
{
	agent_data data(base_data);
	base_data.add();
	data.add();
	base_data.add();
	data.add();
	base_data.add();
	data.add();

	EXPECT_EQ(data.agents_count, 3);

	base_data.remove_at(1);
	data.remove_at(1);
	EXPECT_EQ(data.agents_count, 2);

	base_data.remove_at(0);
	data.remove_at(0);
	EXPECT_EQ(data.agents_count, 1);
}

TEST_F(AgentDataTest, RemoveAtPreservesOtherAgentData)
{
	agent_data data(base_data);

	// Add 3 agents with unique values
	base_data.add();
	data.add();
	data.compartment_types[0] = 7;
	data.cell_ids[0] = 42;
	data.velocities[0] = 0.1;
	data.velocities[1] = 0.2;
	data.velocities[2] = 0.3;

	base_data.add();
	data.add();
	data.compartment_types[1] = 8;
	data.cell_ids[1] = 43;
	data.velocities[3] = 1.1;
	data.velocities[4] = 1.2;
	data.velocities[5] = 1.3;

	base_data.add();
	data.add();
	data.compartment_types[2] = 9;
	data.cell_ids[2] = 44;
	data.velocities[6] = 2.1;
	data.velocities[7] = 2.2;
	data.velocities[8] = 2.3;

	// Remove middle agent
	base_data.remove_at(1);
	data.remove_at(1);

	// Agent 0 should be unchanged
	EXPECT_EQ(data.compartment_types[0], 7);
	EXPECT_EQ(data.cell_ids[0], 42);
	EXPECT_DOUBLE_EQ(data.velocities[0], 0.1);
	EXPECT_DOUBLE_EQ(data.velocities[1], 0.2);
	EXPECT_DOUBLE_EQ(data.velocities[2], 0.3);

	// Agent 2 should now be at index 1 (moved from position 2)
	EXPECT_EQ(data.compartment_types[1], 9);
	EXPECT_EQ(data.cell_ids[1], 44);
	EXPECT_DOUBLE_EQ(data.velocities[3], 2.1);
	EXPECT_DOUBLE_EQ(data.velocities[4], 2.2);
	EXPECT_DOUBLE_EQ(data.velocities[5], 2.3);
}

TEST_F(AgentDataTest, RemoveLastAgentNoMove)
{
	agent_data data(base_data);

	base_data.add();
	data.add();
	data.compartment_types[0] = 11;

	base_data.add();
	data.add();
	data.compartment_types[1] = 22;

	// Remove last agent - no data movement needed
	base_data.remove_at(1);
	data.remove_at(1);

	EXPECT_EQ(data.agents_count, 1);
	EXPECT_EQ(data.compartment_types[0], 11);
}
