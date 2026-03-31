#include <gtest/gtest.h>

#include "micromechanics/cell.h"

using namespace physicore::mechanics::micromechanics;

TEST(CellDataTest, CellDataStructureClear)
{
	cell_data data;
	data.resize(1, 3);
	data.positions[data.cell_offset(0, 0)] = 1.0;
	data.positions[data.cell_offset(0, 1)] = 2.0;
	data.positions[data.cell_offset(0, 2)] = 3.0;
	data.radii[0] = 5.0;
	data.migration_speeds[0] = 2.5;
	data.migration_biases[0] = 0.75;

	data.clear();

	EXPECT_TRUE(data.positions.empty());
	EXPECT_TRUE(data.radii.empty());
	EXPECT_TRUE(data.migration_speeds.empty());
	EXPECT_TRUE(data.migration_biases.empty());
	EXPECT_TRUE(data.agent_counts.empty());
	EXPECT_TRUE(data.pressures.empty());
}

TEST(CellDataTest, PressureAndAgentCount)
{
	cell_data data;
	data.resize(1, 3);
	const cell c0(0, data);

	// Initially zero
	EXPECT_DOUBLE_EQ(c0.pressure(), 0.0);
	EXPECT_EQ(c0.agent_count(), 0);

	// Set pressure directly
	data.pressures[0] = 15.0;
	EXPECT_DOUBLE_EQ(c0.pressure(), 15.0);

	// Set agent count
	data.agent_counts[0] = 5;
	EXPECT_EQ(c0.agent_count(), 5);
}

TEST(CellDataTest, MotilityConfigAccessors)
{
	cell_data data;
	data.resize(1, 3);
	cell c0(0, data);

	c0.migration_speed() = 3.0;
	c0.migration_bias() = 0.4;
	EXPECT_DOUBLE_EQ(data.migration_speeds[0], 3.0);
	EXPECT_DOUBLE_EQ(data.migration_biases[0], 0.4);
}

TEST(CellDataTest, RadiusAccessor)
{
	cell_data data;
	data.resize(1, 3);
	cell c0(0, data);

	c0.radius() = 2.5;
	EXPECT_DOUBLE_EQ(data.radii[0], 2.5);
}

TEST(CellDataTest, ResizePreservesExistingValues)
{
	cell_data data;
	data.resize(1, 3);
	data.radii[0] = 3.0;
	data.migration_speeds[0] = 1.5;

	data.resize(2);
	EXPECT_DOUBLE_EQ(data.radii[0], 3.0);
	EXPECT_DOUBLE_EQ(data.migration_speeds[0], 1.5);
	// New cell has defaults
	EXPECT_DOUBLE_EQ(data.radii[1], 1.0);
	EXPECT_DOUBLE_EQ(data.migration_speeds[1], 0.0);
}
