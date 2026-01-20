#include <gtest/gtest.h>

#include "micromechanics/cell.h"

using namespace physicore::mechanics::micromechanics;
TEST(CellDataTest, CellDataStructureClear)
{
	cell_data data;
	data.resize(1);
	data.positions[data.cell_offset(0, 0)] = 1.0;
	data.positions[data.cell_offset(0, 1)] = 2.0;
	data.positions[data.cell_offset(0, 2)] = 3.0;
	data.volumes[0] = 100.0;
	data.cell_definition_ids[0] = 7;
	data.compartments_count[0] = 12;
	data.migration_speeds[0] = 2.5;
	data.migration_biases[0] = 0.75;
	data.polarizations[data.cell_offset(0, 0)] = 0.1;
	data.polarizations[data.cell_offset(0, 1)] = 0.2;
	data.polarizations[data.cell_offset(0, 2)] = 0.3;
	data.compartment_is_movable[cell_data::compartment_offset(0, 0)] = static_cast<std::uint8_t>(0);

	data.clear();

	EXPECT_TRUE(data.positions.empty());
	EXPECT_TRUE(data.volumes.empty());
	EXPECT_TRUE(data.cell_definition_ids.empty());
	EXPECT_TRUE(data.compartments_count.empty());
	EXPECT_TRUE(data.migration_speeds.empty());
	EXPECT_TRUE(data.migration_biases.empty());
	EXPECT_TRUE(data.polarizations.empty());
	EXPECT_TRUE(data.compartment_is_movable.empty());
}

TEST(CellDataTest, CompartmentPressureMethods)
{
	cell_data data;
	data.resize(1);
	cell c0(0, data);

	// Initially zero
	EXPECT_DOUBLE_EQ(c0.pressure(0), 0.0);

	// Add pressure
	c0.add_pressure(0, 10.0);
	EXPECT_DOUBLE_EQ(c0.pressure(0), 10.0);

	c0.add_pressure(0, 5.0);
	EXPECT_DOUBLE_EQ(c0.pressure(0), 15.0);

	// Different compartment
	c0.add_pressure(1, 20.0);
	EXPECT_DOUBLE_EQ(c0.pressure(1), 20.0);

	// Total cell pressure
	EXPECT_DOUBLE_EQ(c0.total_pressure(), 35.0);
}

TEST(CellDataTest, CompartmentCountMethods)
{
	cell_data data;

	data.resize(1);
	data.compartment_counts[cell_data::compartment_offset(0, 0)] = 3;
	data.compartment_counts[cell_data::compartment_offset(0, 1)] = 2;
	cell c0(0, data);

	EXPECT_EQ(c0.compartment_count(0), 3);
	EXPECT_EQ(c0.compartment_count(1), 2);
	EXPECT_EQ(c0.total_agent_count(), 5);
}

TEST(CellDataTest, MotilityFieldsAccessors)
{
	cell_data data;
	data.resize(1);
	cell c0(0, data);

	c0.migration_speed() = 3.0;
	c0.migration_bias() = 0.4;
	auto pol = c0.polarization();
	ASSERT_EQ(pol.size(), 3u);
	pol[0] = 1.0;
	pol[1] = 0.0;
	pol[2] = 0.0;

	c0.compartment_is_movable(7) = static_cast<std::uint8_t>(0);

	EXPECT_DOUBLE_EQ(data.migration_speeds[0], 3.0);
	EXPECT_DOUBLE_EQ(data.migration_biases[0], 0.4);
	EXPECT_DOUBLE_EQ(data.polarizations[data.cell_offset(0, 0)], 1.0);
	EXPECT_EQ(data.compartment_is_movable[cell_data::compartment_offset(0, 7)], static_cast<std::uint8_t>(0));
}
