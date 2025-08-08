#include <gtest/gtest.h>

#include "base_agent_data.h"

using namespace physicore;

TEST(BaseAgentDataTest, AddAndRemoveAt)
{
	base_agent_data data;
	data.dims = 3;
	data.agents_count = 0;
	data.add();
	EXPECT_EQ(data.agents_count, 1);
	EXPECT_EQ(data.positions.size(), 3);
	data.add();
	EXPECT_EQ(data.agents_count, 2);
	EXPECT_EQ(data.positions.size(), 6);
	data.remove_at(0);
	EXPECT_EQ(data.agents_count, 1);
	EXPECT_EQ(data.positions.size(), 3);
}
