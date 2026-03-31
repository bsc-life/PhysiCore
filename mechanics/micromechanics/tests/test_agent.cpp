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

TEST_F(AgentTest, CompartmentType)
{
	agent test_agent(0, data);
	test_agent.compartment_type() = 7;
	EXPECT_EQ(test_agent.compartment_type(), 7);
}

TEST_F(AgentTest, CellId)
{
	agent test_agent(0, data);
	test_agent.cell_id() = 123;
	EXPECT_EQ(test_agent.cell_id(), 123);
}

TEST_F(AgentTest, Force)
{
	agent test_agent(0, data);
	auto f = test_agent.force();
	ASSERT_EQ(f.size(), 3);
	f[0] = 0.1;
	f[1] = 0.2;
	f[2] = 0.3;
	EXPECT_DOUBLE_EQ(test_agent.force()[0], 0.1);
	EXPECT_DOUBLE_EQ(test_agent.force()[1], 0.2);
	EXPECT_DOUBLE_EQ(test_agent.force()[2], 0.3);
}

TEST_F(AgentTest, SpringAttachments)
{
	agent test_agent(0, data);
	EXPECT_TRUE(test_agent.spring_attachments().empty());
	data.spring_attachments[0].push_back(1);
	data.spring_attachments[0].push_back(2);

	auto att = test_agent.spring_attachments();
	ASSERT_EQ(att.size(), 2);
	EXPECT_EQ(att[0], 1);
	EXPECT_EQ(att[1], 2);
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
	agent0.compartment_type() = 3;
	agent0.cell_id() = 10;

	agent1.velocity()[0] = 2.0;
	agent1.compartment_type() = 4;
	agent1.cell_id() = 11;

	// Verify isolation
	EXPECT_DOUBLE_EQ(agent0.velocity()[0], 1.0);
	EXPECT_EQ(agent0.compartment_type(), 3);
	EXPECT_EQ(agent0.cell_id(), 10);

	EXPECT_DOUBLE_EQ(agent1.velocity()[0], 2.0);
	EXPECT_EQ(agent1.compartment_type(), 4);
	EXPECT_EQ(agent1.cell_id(), 11);
}
