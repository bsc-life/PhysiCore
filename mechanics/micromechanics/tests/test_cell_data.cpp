#include <cmath>
#include <memory>
#include <numbers>

#include <common/base_agent_data.h>
#include <gtest/gtest.h>
#include <openmp_solver/register_solver.h>

#include "micromechanics/agent_container.h"
#include "micromechanics/agent_data.h"
#include "micromechanics/cell_data.h"
#include "micromechanics/environment.h"
#include "micromechanics/solver_registry.h"
#include "micromechanics/uniform_grid_spatial_index.h"

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

// Ensure solver is registered (static library linking doesn't auto-register)
static bool ensure_registered = []() {
	kernels::openmp_solver::attach_to_registry();
	return true;
}();

class CellDataTest : public ::testing::Test
{
protected:
	void SetUp() override
	{
		env = std::make_unique<environment>(0.01);
		auto base_data = std::make_unique<base_agent_data>(3);
		auto mech_data = std::make_unique<agent_data>(*base_data);
		env->agents = std::make_unique<agent_container>(std::move(base_data), std::move(mech_data));
		env->index = std::make_unique<uniform_grid_spatial_index>();

		solver = solver_registry::instance().get("openmp_solver");
	}

	index_t AddAgent(real_t x, real_t y, real_t z, real_t radius, index_t cell_id, std::uint8_t agent_type = 0)
	{
		auto* agent = env->agents->create();
		index_t idx = env->agents->size() - 1;
		agent->position()[0] = x;
		agent->position()[1] = y;
		agent->position()[2] = z;
		agent->radius() = radius;

		auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
		mech_data.cell_ids[idx] = cell_id;
		mech_data.agent_types[idx] = agent_type;
		mech_data.velocities[idx * 3 + 0] = 0.0;
		mech_data.velocities[idx * 3 + 1] = 0.0;
		mech_data.velocities[idx * 3 + 2] = 0.0;

		return idx;
	}

	std::unique_ptr<environment> env;
	solver_ptr solver;
};

TEST_F(CellDataTest, CellDataStructureClear)
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

TEST_F(CellDataTest, CompartmentPressureMethods)
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

TEST_F(CellDataTest, CompartmentCountMethods)
{
	cell_data data;

	EXPECT_EQ(data.get_compartment_count(0, 0), 0);

	data.compartment_counts[{ 0, 0 }] = 3;
	data.compartment_counts[{ 0, 1 }] = 2;

	EXPECT_EQ(data.get_compartment_count(0, 0), 3);
	EXPECT_EQ(data.get_compartment_count(0, 1), 2);
	EXPECT_EQ(data.get_total_agent_count(0), 5);
}

TEST_F(CellDataTest, CalculateCellDataPositions)
{
	// Create a cell with 2 agents at different positions
	AddAgent(10, 0, 0, 5, 0); // cell 0
	AddAgent(20, 0, 0, 5, 0); // cell 0

	solver->initialize(*env);
	solver->calculate_cell_data(*env);

	// Position should be average: (10+20)/2 = 15
	ASSERT_TRUE(env->cells.positions.count(0) > 0);
	EXPECT_DOUBLE_EQ(env->cells.positions[0][0], 15.0);
	EXPECT_DOUBLE_EQ(env->cells.positions[0][1], 0.0);
	EXPECT_DOUBLE_EQ(env->cells.positions[0][2], 0.0);
}

TEST_F(CellDataTest, CalculateCellDataVolumes)
{
	// Create a cell with 2 agents
	AddAgent(0, 0, 0, 10, 0);
	AddAgent(25, 0, 0, 5, 0);

	solver->initialize(*env);
	solver->calculate_cell_data(*env);

	// Volume = sum of 4/3 * pi * r^3
	real_t expected_volume = (4.0 / 3.0) * std::numbers::pi * (10.0 * 10.0 * 10.0 + 5.0 * 5.0 * 5.0);

	ASSERT_TRUE(env->cells.volumes.count(0) > 0);
	EXPECT_NEAR(env->cells.volumes[0], expected_volume, 0.01);
}

TEST_F(CellDataTest, CalculateCellDataVelocities)
{
	index_t idx0 = AddAgent(0, 0, 0, 5, 0);
	index_t idx1 = AddAgent(20, 0, 0, 5, 0);

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
	mech_data.velocities[idx0 * 3 + 0] = 2.0;
	mech_data.velocities[idx1 * 3 + 0] = 4.0;

	solver->initialize(*env);
	solver->calculate_cell_data(*env);

	// Average velocity: (2+4)/2 = 3
	ASSERT_TRUE(env->cells.velocities.count(0) > 0);
	EXPECT_DOUBLE_EQ(env->cells.velocities[0][0], 3.0);

	// Speed = magnitude
	ASSERT_TRUE(env->cells.speeds.count(0) > 0);
	EXPECT_DOUBLE_EQ(env->cells.speeds[0], 3.0);
}

TEST_F(CellDataTest, CalculateCellDataCompartmentCounts)
{
	AddAgent(0, 0, 0, 5, 0, 0);	 // cell 0, type 0
	AddAgent(10, 0, 0, 5, 0, 0); // cell 0, type 0
	AddAgent(20, 0, 0, 5, 0, 1); // cell 0, type 1

	solver->initialize(*env);
	solver->calculate_cell_data(*env);

	EXPECT_EQ(env->cells.get_compartment_count(0, 0), 2);
	EXPECT_EQ(env->cells.get_compartment_count(0, 1), 1);
	EXPECT_EQ(env->cells.get_total_agent_count(0), 3);
}

// ============== Multi-Cell Tests ==============

TEST_F(CellDataTest, MultipleCellsPositions)
{
	// Cell 0 with 2 agents
	AddAgent(0, 0, 0, 5, 0);  // cell 0
	AddAgent(10, 0, 0, 5, 0); // cell 0
	// Cell 1 with 2 agents
	AddAgent(100, 0, 0, 5, 1); // cell 1
	AddAgent(110, 0, 0, 5, 1); // cell 1

	solver->initialize(*env);
	solver->calculate_cell_data(*env);

	// Cell 0 position: average of (0,0,0) and (10,0,0) = (5,0,0)
	ASSERT_TRUE(env->cells.positions.count(0) > 0);
	EXPECT_DOUBLE_EQ(env->cells.positions[0][0], 5.0);

	// Cell 1 position: average of (100,0,0) and (110,0,0) = (105,0,0)
	ASSERT_TRUE(env->cells.positions.count(1) > 0);
	EXPECT_DOUBLE_EQ(env->cells.positions[1][0], 105.0);
}

TEST_F(CellDataTest, MultipleCellsVolumes)
{
	// Cell 0 with 1 agent (radius 10)
	AddAgent(0, 0, 0, 10, 0);
	// Cell 1 with 2 agents (radius 5 each)
	AddAgent(100, 0, 0, 5, 1);
	AddAgent(110, 0, 0, 5, 1);

	solver->initialize(*env);
	solver->calculate_cell_data(*env);

	real_t vol_r10 = (4.0 / 3.0) * std::numbers::pi * 10.0 * 10.0 * 10.0;
	real_t vol_r5 = (4.0 / 3.0) * std::numbers::pi * 5.0 * 5.0 * 5.0;

	ASSERT_TRUE(env->cells.volumes.count(0) > 0);
	EXPECT_NEAR(env->cells.volumes[0], vol_r10, 0.01);

	ASSERT_TRUE(env->cells.volumes.count(1) > 0);
	EXPECT_NEAR(env->cells.volumes[1], 2 * vol_r5, 0.01);
}

TEST_F(CellDataTest, MultipleCellsCompartmentCounts)
{
	// Cell 0: 2 agents of type 0, 1 agent of type 1
	AddAgent(0, 0, 0, 5, 0, 0);
	AddAgent(10, 0, 0, 5, 0, 0);
	AddAgent(20, 0, 0, 5, 0, 1);
	// Cell 1: 1 agent of type 0, 2 agents of type 1
	AddAgent(100, 0, 0, 5, 1, 0);
	AddAgent(110, 0, 0, 5, 1, 1);
	AddAgent(120, 0, 0, 5, 1, 1);

	solver->initialize(*env);
	solver->calculate_cell_data(*env);

	// Cell 0 compartments
	EXPECT_EQ(env->cells.get_compartment_count(0, 0), 2);
	EXPECT_EQ(env->cells.get_compartment_count(0, 1), 1);
	EXPECT_EQ(env->cells.get_total_agent_count(0), 3);

	// Cell 1 compartments
	EXPECT_EQ(env->cells.get_compartment_count(1, 0), 1);
	EXPECT_EQ(env->cells.get_compartment_count(1, 1), 2);
	EXPECT_EQ(env->cells.get_total_agent_count(1), 3);
}

TEST_F(CellDataTest, MultipleCellsVelocities)
{
	// Cell 0: 2 agents with different velocities
	index_t idx0 = AddAgent(0, 0, 0, 5, 0);
	index_t idx1 = AddAgent(10, 0, 0, 5, 0);
	// Cell 1: 1 agent
	index_t idx2 = AddAgent(100, 0, 0, 5, 1);

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
	mech_data.velocities[idx0 * 3 + 0] = 4.0;
	mech_data.velocities[idx1 * 3 + 0] = 6.0;
	mech_data.velocities[idx2 * 3 + 0] = 10.0;

	solver->initialize(*env);
	solver->calculate_cell_data(*env);

	// Cell 0: average velocity = (4+6)/2 = 5
	ASSERT_TRUE(env->cells.velocities.count(0) > 0);
	EXPECT_DOUBLE_EQ(env->cells.velocities[0][0], 5.0);
	EXPECT_DOUBLE_EQ(env->cells.speeds[0], 5.0);

	// Cell 1: velocity = 10
	ASSERT_TRUE(env->cells.velocities.count(1) > 0);
	EXPECT_DOUBLE_EQ(env->cells.velocities[1][0], 10.0);
	EXPECT_DOUBLE_EQ(env->cells.speeds[1], 10.0);
}
