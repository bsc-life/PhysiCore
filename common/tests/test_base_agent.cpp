#include <gtest/gtest.h>

#include "base_agent.h"
#include "base_agent_data.h"

using namespace physicore;

TEST(BaseAgentTest, GetPosition)
{
	base_agent_data data;
	data.dims = 2;
	data.agents_count = 0;
	data.add();
	data.positions[0] = 1.0;
	data.positions[1] = 2.0;
	base_agent agent(0, data);
	auto pos = agent.position();
	ASSERT_EQ(pos.size(), 2);
	EXPECT_EQ(pos[0], 1.0);
	EXPECT_EQ(pos[1], 2.0);
}

TEST(BaseAgentTest, GetPosition3D)
{
	base_agent_data data;
	data.dims = 3;
	data.agents_count = 0;
	data.add();
	data.positions[0] = 3.0;
	data.positions[1] = 4.0;
	data.positions[2] = 5.0;
	base_agent agent(0, data);
	auto pos = agent.position();
	ASSERT_EQ(pos.size(), 3);
	EXPECT_EQ(pos[0], 3.0);
	EXPECT_EQ(pos[1], 4.0);
	EXPECT_EQ(pos[2], 5.0);
}

TEST(BaseAgentTest, MultipleAgentsGetPosition)
{
	base_agent_data data;
	data.dims = 2;
	data.agents_count = 0;
	data.add(); // agent 0
	data.add(); // agent 1
	data.positions[0] = 10.0;
	data.positions[1] = 20.0;
	data.positions[2] = 30.0;
	data.positions[3] = 40.0;
	base_agent agent0(0, data);
	base_agent agent1(1, data);
	auto pos0 = agent0.position();
	auto pos1 = agent1.position();
	ASSERT_EQ(pos0.size(), 2);
	ASSERT_EQ(pos1.size(), 2);
	EXPECT_EQ(pos0[0], 10.0);
	EXPECT_EQ(pos0[1], 20.0);
	EXPECT_EQ(pos1[0], 30.0);
	EXPECT_EQ(pos1[1], 40.0);
}
