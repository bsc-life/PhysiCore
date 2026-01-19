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

class EnvironmentTest : public ::testing::Test
{
protected:
	static std::unique_ptr<environment> create_test_environment()
	{
		auto env = std::make_unique<environment>(0.01);
		auto base_data = std::make_unique<base_agent_data>(3);
		auto mech_data = std::make_unique<agent_data>(*base_data);
		env->agents = std::make_unique<agent_container>(std::move(base_data), std::move(mech_data));
		env->index = std::make_unique<uniform_grid_spatial_index>();
		return env;
	}
};

TEST_F(EnvironmentTest, RunSingleTimestep)
{
	auto env = create_test_environment();
	env->run_single_timestep();
	SUCCEED();
}

TEST_F(EnvironmentTest, RunMultipleTimestepsWithForces)
{
	auto env = create_test_environment();

	// Set up the solver
	env->solver_ = solver_registry::instance().get("openmp_solver");
	ASSERT_NE(env->solver_, nullptr);
	env->solver_->initialize(*env);

	// Add two overlapping agents
	auto* agent0 = env->agents->create();
	agent0->position()[0] = 0.0;
	agent0->radius() = 10.0;
	agent0->is_movable() = 1;
	agent0->cell_cell_repulsion_strength() = 10.0;
	agent0->relative_maximum_adhesion_distance() = 1.5;

	auto* agent1 = env->agents->create();
	agent1->position()[0] = 15.0;
	agent1->radius() = 10.0;
	agent1->is_movable() = 1;
	agent1->cell_cell_repulsion_strength() = 10.0;
	agent1->relative_maximum_adhesion_distance() = 1.5;

	real_t const initial_x0 = agent0->position()[0];
	real_t const initial_x1 = agent1->position()[0];

	// Run several timesteps
	for (int i = 0; i < 10; ++i)
	{
		env->run_single_timestep();
	}

	// Agents should have moved apart due to repulsion
	EXPECT_LT(agent0->position()[0], initial_x0);
	EXPECT_GT(agent1->position()[0], initial_x1);
}
