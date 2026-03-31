#include <memory>

#include <common/base_agent_data.h>
#include <gtest/gtest.h>

#include "micromechanics/agent_container.h"
#include "micromechanics/agent_data.h"
#include "micromechanics/environment.h"
#include "micromechanics/solver.h"
#include "micromechanics/solver_registry.h"
#include "micromechanics/uniform_grid_spatial_index.h"

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

namespace {

/// Mock solver that tracks which methods were called.
class tracking_solver : public solver
{
public:
	bool initialized = false;
	bool neighbors_updated = false;
	bool forces_updated = false;
	bool cell_data_calculated = false;
	bool motility_updated = false;
	bool membrane_updated = false;
	bool springs_updated = false;
	bool positions_updated = false;

	void initialize(environment& /*e*/) override { initialized = true; }
	void update_cell_neighbors(environment& /*e*/) override { neighbors_updated = true; }
	void update_cell_forces(environment& /*e*/) override { forces_updated = true; }
	void calculate_cell_data(environment& /*e*/) override { cell_data_calculated = true; }
	void update_motility(environment& /*e*/) override { motility_updated = true; }
	void update_basement_membrane_interactions(environment& /*e*/) override { membrane_updated = true; }
	void update_spring_attachments(environment& /*e*/) override { springs_updated = true; }
	void update_positions(environment& /*e*/) override { positions_updated = true; }
};

} // namespace

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

TEST_F(EnvironmentTest, ConstructorDefaults)
{
	const environment env(0.05);

	EXPECT_DOUBLE_EQ(env.timestep, 0.05);
	EXPECT_NE(env.agents, nullptr);
	EXPECT_NE(env.index, nullptr);
	EXPECT_EQ(env.solver_, nullptr);

	// Domain defaults
	EXPECT_DOUBLE_EQ(env.domain_min[0], -500.0);
	EXPECT_DOUBLE_EQ(env.domain_min[1], -500.0);
	EXPECT_DOUBLE_EQ(env.domain_min[2], -500.0);
	EXPECT_DOUBLE_EQ(env.domain_max[0], 500.0);
	EXPECT_DOUBLE_EQ(env.domain_max[1], 500.0);
	EXPECT_DOUBLE_EQ(env.domain_max[2], 500.0);
}

TEST_F(EnvironmentTest, RunSingleTimestepWithoutSolver)
{
	auto env = create_test_environment();
	// solver_ is null — should be a no-op, not crash
	env->run_single_timestep();
	SUCCEED();
}

TEST_F(EnvironmentTest, InitializeSolverFromRegistry)
{
	// Register our tracking solver under a test-specific name
	registry_adder<tracking_solver> const adder("test_tracking_solver_init");

	auto env = create_test_environment();
	env->params.solver_name = "test_tracking_solver_init";
	env->initialize_solver();

	ASSERT_NE(env->solver_, nullptr);
	// Verify initialize was called
	auto* ts = dynamic_cast<tracking_solver*>(env->solver_.get());
	ASSERT_NE(ts, nullptr);
	EXPECT_TRUE(ts->initialized);
}

TEST_F(EnvironmentTest, InitializeSolverUnknownName)
{
	auto env = create_test_environment();
	env->params.solver_name = "nonexistent_solver_xyz";

#ifdef NDEBUG
	env->initialize_solver();
	EXPECT_EQ(env->solver_, nullptr);
#else
	SUCCEED(); // In debug mode, get() may assert — skip
#endif
}

TEST_F(EnvironmentTest, RunTimestepMotilityOnly)
{
	auto env = create_test_environment();
	auto ts = std::make_unique<tracking_solver>();
	auto* ts_ptr = ts.get();
	env->solver_ = std::move(ts);

	env->params.enable_motility = true;
	env->params.enable_basement_membrane = false;
	env->params.enable_spring_attachments = false;

	env->run_single_timestep();

	EXPECT_TRUE(ts_ptr->neighbors_updated);
	EXPECT_TRUE(ts_ptr->forces_updated);
	EXPECT_TRUE(ts_ptr->cell_data_calculated);
	EXPECT_TRUE(ts_ptr->motility_updated);
	EXPECT_FALSE(ts_ptr->membrane_updated);
	EXPECT_FALSE(ts_ptr->springs_updated);
	EXPECT_TRUE(ts_ptr->positions_updated);
}

TEST_F(EnvironmentTest, RunTimestepBasementMembraneOnly)
{
	auto env = create_test_environment();
	auto ts = std::make_unique<tracking_solver>();
	auto* ts_ptr = ts.get();
	env->solver_ = std::move(ts);

	env->params.enable_motility = false;
	env->params.enable_basement_membrane = true;
	env->params.enable_spring_attachments = false;

	env->run_single_timestep();

	EXPECT_FALSE(ts_ptr->motility_updated);
	EXPECT_TRUE(ts_ptr->membrane_updated);
	EXPECT_FALSE(ts_ptr->springs_updated);
}

TEST_F(EnvironmentTest, RunTimestepSpringAttachmentsOnly)
{
	auto env = create_test_environment();
	auto ts = std::make_unique<tracking_solver>();
	auto* ts_ptr = ts.get();
	env->solver_ = std::move(ts);

	env->params.enable_motility = false;
	env->params.enable_basement_membrane = false;
	env->params.enable_spring_attachments = true;

	env->run_single_timestep();

	EXPECT_FALSE(ts_ptr->motility_updated);
	EXPECT_FALSE(ts_ptr->membrane_updated);
	EXPECT_TRUE(ts_ptr->springs_updated);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_F(EnvironmentTest, RunTimestepAllFeaturesEnabled)
{
	auto env = create_test_environment();
	auto ts = std::make_unique<tracking_solver>();
	auto* ts_ptr = ts.get();
	env->solver_ = std::move(ts);

	env->params.enable_motility = true;
	env->params.enable_basement_membrane = true;
	env->params.enable_spring_attachments = true;

	env->run_single_timestep();

	EXPECT_TRUE(ts_ptr->neighbors_updated);
	EXPECT_TRUE(ts_ptr->forces_updated);
	EXPECT_TRUE(ts_ptr->cell_data_calculated);
	EXPECT_TRUE(ts_ptr->motility_updated);
	EXPECT_TRUE(ts_ptr->membrane_updated);
	EXPECT_TRUE(ts_ptr->springs_updated);
	EXPECT_TRUE(ts_ptr->positions_updated);
}

TEST_F(EnvironmentTest, RunTimestepAllFeaturesDisabled)
{
	auto env = create_test_environment();
	auto ts = std::make_unique<tracking_solver>();
	auto* ts_ptr = ts.get();
	env->solver_ = std::move(ts);

	env->params.enable_motility = false;
	env->params.enable_basement_membrane = false;
	env->params.enable_spring_attachments = false;

	env->run_single_timestep();

	// Core pipeline always runs
	EXPECT_TRUE(ts_ptr->neighbors_updated);
	EXPECT_TRUE(ts_ptr->forces_updated);
	EXPECT_TRUE(ts_ptr->cell_data_calculated);
	EXPECT_TRUE(ts_ptr->positions_updated);

	// Optional features skipped
	EXPECT_FALSE(ts_ptr->motility_updated);
	EXPECT_FALSE(ts_ptr->membrane_updated);
	EXPECT_FALSE(ts_ptr->springs_updated);
}

TEST_F(EnvironmentTest, SerializeStateDoesNotCrash)
{
	auto env = create_test_environment();
	env->serialize_state(1.0);
	env->serialize_state(0.0);
	env->serialize_state(-1.0);
	SUCCEED();
}
