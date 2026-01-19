#include <gtest/gtest.h>
#include "micromechanics/cell_data.h"
using namespace physicore::mechanics::micromechanics;
TEST(CellDataTest, CellDataStructureClear)
{
	cell_data data;
	data.positions[0] = { 1.0, 2.0, 3.0 };
	data.volumes[0] = 100.0;
	data.speeds[0] = 5.0;

	data.clear();

	EXPECT_TRUE(data.positions.empty());
	EXPECT_TRUE(data.volumes.empty());
	EXPECT_TRUE(data.speeds.empty());
}

TEST(CellDataTest, CompartmentPressureMethods)
{
	cell_data data;

	// Initially zero
	EXPECT_DOUBLE_EQ(data.get_pressure(0, 0), 0.0);

	// Add pressure
	data.add_pressure(0, 0, 10.0);
	EXPECT_DOUBLE_EQ(data.get_pressure(0, 0), 10.0);

	data.add_pressure(0, 0, 5.0);
	EXPECT_DOUBLE_EQ(data.get_pressure(0, 0), 15.0);

	// Different compartment
	data.add_pressure(0, 1, 20.0);
	EXPECT_DOUBLE_EQ(data.get_pressure(0, 1), 20.0);

	// Total cell pressure
	EXPECT_DOUBLE_EQ(data.get_total_cell_pressure(0), 35.0);
}

TEST(CellDataTest, CompartmentCountMethods)
{
	cell_data data;

	EXPECT_EQ(data.get_compartment_count(0, 0), 0);

	data.compartment_counts[{ 0, 0 }] = 3;
	data.compartment_counts[{ 0, 1 }] = 2;

	EXPECT_EQ(data.get_compartment_count(0, 0), 3);
	EXPECT_EQ(data.get_compartment_count(0, 1), 2);
	EXPECT_EQ(data.get_total_agent_count(0), 5);
}
