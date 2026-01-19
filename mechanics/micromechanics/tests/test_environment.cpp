#include <filesystem>
#include <fstream>
#include <memory>
#include <string>

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
	void SetUp() override
	{
		// Create a temporary output directory for tests
		test_output_dir = std::filesystem::temp_directory_path() / "vtk_mechanics_test";
		std::filesystem::create_directories(test_output_dir);

		// Clean up any existing files
		if (std::filesystem::exists(test_output_dir))
		{
			std::filesystem::remove_all(test_output_dir);
		}
		std::filesystem::create_directories(test_output_dir);
	}

	void TearDown() override
	{
		// Clean up test files
		if (std::filesystem::exists(test_output_dir))
		{
			std::filesystem::remove_all(test_output_dir);
		}
	}

	static std::unique_ptr<environment> create_test_environment()
	{
		auto env = std::make_unique<environment>(0.01);
		auto base_data = std::make_unique<base_agent_data>(3);
		auto mech_data = std::make_unique<agent_data>(*base_data);
		env->agents = std::make_unique<agent_container>(std::move(base_data), std::move(mech_data));
		env->index = std::make_unique<uniform_grid_spatial_index>();
		return env;
	}

	std::filesystem::path test_output_dir;
};

TEST_F(EnvironmentTest, RunSingleTimestep)
{
	auto env = create_test_environment();
	env->run_single_timestep();
	SUCCEED();
}

TEST_F(EnvironmentTest, SerializeCreatesFiles)
{
	auto env = create_test_environment();

	// Add agents to serialize
	for (int i = 0; i < 5; ++i)
	{
		auto* agent = env->agents->create();
		agent->position()[0] = i * 10.0;
		agent->position()[1] = 0.0;
		agent->position()[2] = 0.0;
		agent->radius() = 5.0;
	}

	env->serialize_state(0.0);

	// Check that VTK file is created
	auto vtk_dir = std::filesystem::path("output") / "vtk_mechanics";
	auto vtu_file = vtk_dir / "mechanics_000000.vtu";
	EXPECT_TRUE(std::filesystem::exists(vtu_file));

	// Check that PVD file is created
	auto pvd_file = std::filesystem::path("output") / "mechanics.pvd";
	EXPECT_TRUE(std::filesystem::exists(pvd_file));
}

TEST_F(EnvironmentTest, SerializeMultipleTimes)
{
	auto env = create_test_environment();

	// Add an agent
	auto* agent = env->agents->create();
	agent->position()[0] = 0.0;
	agent->radius() = 5.0;

	// Serialize multiple times
	for (int i = 0; i < 3; ++i)
	{
		env->serialize_state(i * 0.1);
	}

	auto vtk_dir = std::filesystem::path("output") / "vtk_mechanics";

	// Check that multiple VTK files are created
	EXPECT_TRUE(std::filesystem::exists(vtk_dir / "mechanics_000000.vtu"));
	EXPECT_TRUE(std::filesystem::exists(vtk_dir / "mechanics_000001.vtu"));
	EXPECT_TRUE(std::filesystem::exists(vtk_dir / "mechanics_000002.vtu"));
}

TEST_F(EnvironmentTest, PvdFileContainsCorrectEntries)
{
	auto env = create_test_environment();

	auto* agent = env->agents->create();
	agent->position()[0] = 0.0;
	agent->radius() = 5.0;

	// Serialize twice with different times
	env->serialize_state(0.1);
	env->serialize_state(0.2);

	// Read PVD file content
	auto pvd_file = std::filesystem::path("output") / "mechanics.pvd";
	ASSERT_TRUE(std::filesystem::exists(pvd_file));

	std::ifstream file(pvd_file);
	std::string const content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

	// Check XML structure
	EXPECT_TRUE(content.find("<?xml version=\"1.0\"?>") != std::string::npos);
	EXPECT_TRUE(content.find("<VTKFile type=\"Collection\"") != std::string::npos);
	EXPECT_TRUE(content.find("<Collection>") != std::string::npos);
	EXPECT_TRUE(content.find("</Collection>") != std::string::npos);
	EXPECT_TRUE(content.find("</VTKFile>") != std::string::npos);

	// Check timestep entries
	EXPECT_TRUE(content.find("timestep=\"0.1") != std::string::npos);
	EXPECT_TRUE(content.find("timestep=\"0.2") != std::string::npos);

	// Check file references
	EXPECT_TRUE(content.find("mechanics_000000.vtu") != std::string::npos);
	EXPECT_TRUE(content.find("mechanics_000001.vtu") != std::string::npos);
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
