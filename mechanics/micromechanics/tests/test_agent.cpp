#include <common/base_agent_data.h>
#include <gtest/gtest.h>

#include "micromechanics/agent.h"
#include "micromechanics/agent_data.h"

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

class AgentTest : public ::testing::Test
{
protected:
	void SetUp() override
	{
		base_data.dims = 3;
		base_data.agents_count = 0;

		// Add first agent
		base_data.add();
		data.add();
	}

	base_agent_data base_data;
	agent_data data = agent_data(base_data);
};

TEST_F(AgentTest, Velocity)
{
	agent test_agent(0, data);
	auto vel = test_agent.velocity();
	ASSERT_EQ(vel.size(), 3);

	vel[0] = 1.0;
	vel[1] = 2.0;
	vel[2] = 3.0;
	EXPECT_DOUBLE_EQ(test_agent.velocity()[0], 1.0);
	EXPECT_DOUBLE_EQ(test_agent.velocity()[1], 2.0);
	EXPECT_DOUBLE_EQ(test_agent.velocity()[2], 3.0);
}

TEST_F(AgentTest, PreviousVelocity)
{
	agent test_agent(0, data);
	auto prev_vel = test_agent.previous_velocity();
	ASSERT_EQ(prev_vel.size(), 3);

	prev_vel[0] = 0.5;
	prev_vel[1] = 1.5;
	prev_vel[2] = 2.5;
	EXPECT_DOUBLE_EQ(test_agent.previous_velocity()[0], 0.5);
	EXPECT_DOUBLE_EQ(test_agent.previous_velocity()[1], 1.5);
	EXPECT_DOUBLE_EQ(test_agent.previous_velocity()[2], 2.5);
}

TEST_F(AgentTest, Position)
{
	agent test_agent(0, data);
	auto pos = test_agent.position();
	ASSERT_EQ(pos.size(), 3);

	pos[0] = 10.0;
	pos[1] = 20.0;
	pos[2] = 30.0;
	EXPECT_DOUBLE_EQ(test_agent.position()[0], 10.0);
	EXPECT_DOUBLE_EQ(test_agent.position()[1], 20.0);
	EXPECT_DOUBLE_EQ(test_agent.position()[2], 30.0);
}

TEST_F(AgentTest, Radius)
{
	agent test_agent(0, data);
	test_agent.radius() = 8.5;
	EXPECT_DOUBLE_EQ(test_agent.radius(), 8.5);
}

TEST_F(AgentTest, IsMovable)
{
	agent test_agent(0, data);
	// Default should be 1 (movable)
	EXPECT_EQ(test_agent.is_movable(), 1);

	test_agent.is_movable() = 0;
	EXPECT_EQ(test_agent.is_movable(), 0);
}

TEST_F(AgentTest, CellCellAdhesionStrength)
{
	agent test_agent(0, data);
	test_agent.cell_cell_adhesion_strength() = 0.4;
	EXPECT_DOUBLE_EQ(test_agent.cell_cell_adhesion_strength(), 0.4);
}

TEST_F(AgentTest, CellCellRepulsionStrength)
{
	agent test_agent(0, data);
	test_agent.cell_cell_repulsion_strength() = 10.0;
	EXPECT_DOUBLE_EQ(test_agent.cell_cell_repulsion_strength(), 10.0);
}

TEST_F(AgentTest, RelativeMaximumAdhesionDistance)
{
	agent test_agent(0, data);
	test_agent.relative_maximum_adhesion_distance() = 1.25;
	EXPECT_DOUBLE_EQ(test_agent.relative_maximum_adhesion_distance(), 1.25);
}

TEST_F(AgentTest, Neighbors)
{
	agent test_agent(0, data);
	// Initially empty
	EXPECT_TRUE(test_agent.neighbors().empty());

	// Add neighbors to the underlying data
	data.neighbors[0].push_back(1);
	data.neighbors[0].push_back(2);

	auto neighbors = test_agent.neighbors();
	ASSERT_EQ(neighbors.size(), 2);
	EXPECT_EQ(neighbors[0], 1);
	EXPECT_EQ(neighbors[1], 2);
}

TEST_F(AgentTest, MultipleAgents)
{
	// Add second agent
	base_data.add();
	data.add();

	agent agent0(0, data);
	agent agent1(1, data);

	// Set unique values
	agent0.velocity()[0] = 1.0;
	agent0.radius() = 5.0;
	agent0.cell_cell_repulsion_strength() = 10.0;

	agent1.velocity()[0] = 2.0;
	agent1.radius() = 7.0;
	agent1.cell_cell_repulsion_strength() = 15.0;

	// Verify isolation
	EXPECT_DOUBLE_EQ(agent0.velocity()[0], 1.0);
	EXPECT_DOUBLE_EQ(agent0.radius(), 5.0);
	EXPECT_DOUBLE_EQ(agent0.cell_cell_repulsion_strength(), 10.0);

	EXPECT_DOUBLE_EQ(agent1.velocity()[0], 2.0);
	EXPECT_DOUBLE_EQ(agent1.radius(), 7.0);
	EXPECT_DOUBLE_EQ(agent1.cell_cell_repulsion_strength(), 15.0);
}

