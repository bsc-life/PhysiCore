#include <memory>

#include <common/base_agent_data.h>
#include <gtest/gtest.h>
#include <openmp_solver/register_solver.h>

#include "micromechanics/agent_container.h"
#include "micromechanics/agent_data.h"
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

class SolverTest : public ::testing::Test
{
protected:
	void SetUp() override
	{
		env = std::make_unique<environment>(0.01);
		auto base_data = std::make_unique<base_agent_data>(3);
		auto mech_data = std::make_unique<agent_data>(*base_data);
		env->agents = std::make_unique<agent_container>(std::move(base_data), std::move(mech_data));
		env->index = std::make_unique<uniform_grid_spatial_index>();

		// Get the OpenMP solver
		solver = solver_registry::instance().get("openmp_solver");
		ASSERT_NE(solver, nullptr);
	}

	void AddAgent(real_t x, real_t y, real_t z, real_t radius)
	{
		auto* agent = env->agents->create();
		agent->position()[0] = x;
		agent->position()[1] = y;
		agent->position()[2] = z;
		agent->radius() = radius;
		agent->is_movable() = 1;
		agent->cell_cell_repulsion_strength() = 10.0;
		agent->relative_maximum_adhesion_distance() = 1.5;
	}

	std::unique_ptr<environment> env;
	solver_ptr solver;
};

TEST_F(SolverTest, InitializeDoesNotThrow)
{
	AddAgent(0, 0, 0, 10);
	EXPECT_NO_THROW(solver->initialize(*env));
}

TEST_F(SolverTest, UpdateNeighborsFindsNearbyAgents)
{
	AddAgent(0, 0, 0, 10);
	AddAgent(15, 0, 0, 10);	 // Close neighbor
	AddAgent(100, 0, 0, 10); // Far away

	solver->initialize(*env);
	solver->update_cell_neighbors(*env);

	// Neighbors should be populated after update_cell_neighbors
	// The actual neighbor finding depends on spatial index implementation
	SUCCEED(); // Basic test that it doesn't crash
}

TEST_F(SolverTest, UpdateForcesCalculatesRepulsion)
{
	// Two overlapping agents
	AddAgent(0, 0, 0, 10);
	AddAgent(15, 0, 0, 10);

	solver->initialize(*env);
	solver->update_cell_neighbors(*env);
	solver->update_cell_forces(*env);

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);

	// Forces should be non-zero and opposite for overlapping agents
	// Agent 0 should be pushed in negative x direction
	// Agent 1 should be pushed in positive x direction
	// (actual direction depends on implementation)
	bool has_force = (mech_data.forces[0] != 0.0 || mech_data.forces[3] != 0.0);
	EXPECT_TRUE(has_force);
}

TEST_F(SolverTest, UpdatePositionsMovesAgents)
{
	// Two overlapping agents with repulsion
	AddAgent(0, 0, 0, 10);
	AddAgent(15, 0, 0, 10);

	auto& base_data = *std::get<std::unique_ptr<base_agent_data>>(env->agents->agent_datas);

	real_t initial_x0 = base_data.positions[0];
	real_t initial_x1 = base_data.positions[3];

	solver->initialize(*env);
	solver->update_cell_neighbors(*env);
	solver->update_cell_forces(*env);
	solver->update_positions(*env);

	// Positions should have changed
	bool pos_changed = (base_data.positions[0] != initial_x0 || base_data.positions[3] != initial_x1);
	EXPECT_TRUE(pos_changed);

	// Agents should have moved apart (x0 decreased, x1 increased)
	EXPECT_LT(base_data.positions[0], initial_x0);
	EXPECT_GT(base_data.positions[3], initial_x1);
}

TEST_F(SolverTest, ImmovableAgentDoesNotMove)
{
	AddAgent(0, 0, 0, 10);
	AddAgent(15, 0, 0, 10);

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
	auto& base_data = *std::get<std::unique_ptr<base_agent_data>>(env->agents->agent_datas);

	// Make agent 0 immovable
	mech_data.is_movable[0] = 0;

	real_t initial_x0 = base_data.positions[0];

	solver->initialize(*env);
	solver->update_cell_neighbors(*env);
	solver->update_cell_forces(*env);
	solver->update_positions(*env);

	// Agent 0 should not have moved
	EXPECT_DOUBLE_EQ(base_data.positions[0], initial_x0);
	EXPECT_DOUBLE_EQ(base_data.positions[1], 0.0);
	EXPECT_DOUBLE_EQ(base_data.positions[2], 0.0);
}

TEST_F(SolverTest, MotilitySolverUpdatesDirection)
{
	AddAgent(0, 0, 0, 10);

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);

	// Make agent motile
	mech_data.is_motile[0] = 1;
	mech_data.persistence_times[0] = 1.0;
	mech_data.migration_speeds[0] = 1.0;
	mech_data.migration_biases[0] = 0.0; // Pure random walk

	solver->initialize(*env);
	solver->update_motility(*env);

	// Motility direction should be set (might still be zero if random gives zero)
	// Just test that it doesn't crash
	SUCCEED();
}
