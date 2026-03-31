#include <gtest/gtest.h>

#include "micromechanics/cell_aggregation.h"

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

namespace {
void add_agent(base_agent_data& base, agent_data& mech)
{
	base.add();
	mech.add();
}
} // namespace

TEST(CellAggregationTest, AggregatesCOMCountsAndPressureProxy)
{
	base_agent_data base(3);
	agent_data mech(base);

	cell_data cells;
	cells.resize(2, 3);

	// Pre-fill to ensure aggregation resets
	cells.positions[cells.cell_offset(0, 0)] = 123.0;
	cells.pressures[0] = 99.0;

	add_agent(base, mech); // a0
	add_agent(base, mech); // a1
	add_agent(base, mech); // a2
	add_agent(base, mech); // a3 (ignored)

	// a0: cell 0
	mech.cell_ids[0] = 0;
	mech.compartment_types[0] = 0;
	base.positions[0 * 3 + 0] = 1.0;
	base.positions[0 * 3 + 1] = 0.0;
	base.positions[0 * 3 + 2] = 0.0;
	mech.velocities[0 * 3 + 0] = 1.0;
	mech.velocities[0 * 3 + 1] = 0.0;
	mech.velocities[0 * 3 + 2] = 0.0;
	mech.forces[0 * 3 + 0] = 3.0;
	mech.forces[0 * 3 + 1] = 4.0;
	mech.forces[0 * 3 + 2] = 0.0; // ||F|| = 5

	// a1: cell 0
	mech.cell_ids[1] = 0;
	mech.compartment_types[1] = 1;
	base.positions[1 * 3 + 0] = 3.0;
	base.positions[1 * 3 + 1] = 0.0;
	base.positions[1 * 3 + 2] = 0.0;
	mech.velocities[1 * 3 + 0] = 1.0;
	mech.velocities[1 * 3 + 1] = 0.0;
	mech.velocities[1 * 3 + 2] = 0.0;
	mech.forces[1 * 3 + 0] = 0.0;
	mech.forces[1 * 3 + 1] = 0.0;
	mech.forces[1 * 3 + 2] = 0.0;

	// a2: cell 1
	mech.cell_ids[2] = 1;
	mech.compartment_types[2] = 0;
	base.positions[2 * 3 + 0] = 0.0;
	base.positions[2 * 3 + 1] = 2.0;
	base.positions[2 * 3 + 2] = 0.0;
	mech.velocities[2 * 3 + 0] = 0.0;
	mech.velocities[2 * 3 + 1] = 2.0;
	mech.velocities[2 * 3 + 2] = 0.0;
	mech.forces[2 * 3 + 0] = 0.0;
	mech.forces[2 * 3 + 1] = 0.0;
	mech.forces[2 * 3 + 2] = 6.0; // ||F|| = 6

	// a3: ignored (invalid cell_id)
	mech.cell_ids[3] = -1;
	mech.compartment_types[3] = 0;

	aggregate_cell_data_from_agents(base, mech, cells);

	// cell 0 COM position: ( (1,0,0) + (3,0,0) ) / 2 = (2,0,0)
	EXPECT_DOUBLE_EQ(cells.positions[cells.cell_offset(0, 0)], 2.0);
	EXPECT_DOUBLE_EQ(cells.positions[cells.cell_offset(0, 1)], 0.0);
	EXPECT_DOUBLE_EQ(cells.positions[cells.cell_offset(0, 2)], 0.0);

	// cell 0 COM velocity: (1,0,0)
	EXPECT_DOUBLE_EQ(cells.velocities[cells.cell_offset(0, 0)], 1.0);

	// cell 0 has 2 agents, cell 1 has 1
	EXPECT_EQ(cells.agent_counts[0], 2);
	EXPECT_EQ(cells.agent_counts[1], 1);

	// cell 0 pressure: ||F_a0|| + ||F_a1|| = 5 + 0 = 5
	EXPECT_NEAR(cells.pressures[0], 5.0, 1e-12);

	// cell 1 position & pressure
	EXPECT_DOUBLE_EQ(cells.positions[cells.cell_offset(1, 0)], 0.0);
	EXPECT_DOUBLE_EQ(cells.positions[cells.cell_offset(1, 1)], 2.0);
	EXPECT_NEAR(cells.pressures[1], 6.0, 1e-12);
}
