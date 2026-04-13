#include <memory>
#include <stdexcept>

#include <gtest/gtest.h>

#include "physicell/environment.h"
#include "physicell/serializer.h"
#include "physicell/solver.h"

namespace physicore::mechanics::physicell::tests {

using namespace physicore;
using namespace physicore::mechanics::physicell;

// Mock solver for testing
class MockSolver : public solver
{
public:
	void initialize(environment& /* e */) override {}
	void solve(environment& /* e */, index_t /* iterations */) override { solve_called = true; }
	bool solve_called = false;
};

// Mock serializer for testing
class MockSerializer : public serializer
{
public:
	void serialize(real_t /* current_time */) override { serialize_called = true; }
	bool serialize_called = false;
};

class EnvironmentTest : public ::testing::Test
{
};

// Test constructor and basic initialization
TEST_F(EnvironmentTest, ConstructorInitializesEnvironment)
{
	environment env(0.1, 3, 2, 3);
	EXPECT_DOUBLE_EQ(env.mechanics_timestep, 0.1);
	EXPECT_FALSE(env.automated_spring_adhesion);
	EXPECT_TRUE(env.virtual_wall_at_domain_edges);
	ASSERT_NE(env.agents, nullptr);
	EXPECT_EQ(env.agents->size(), 0);
}

// Test run_single_timestep calls solver and serialize_state calls serializer
TEST_F(EnvironmentTest, RunSingleTimestepAndSerializeState)
{
	environment env(0.1, 3, 2, 3);

	auto mock_solver = std::make_unique<MockSolver>();
	auto mock_serializer = std::make_unique<MockSerializer>();

	auto* solver_ptr = mock_solver.get();
	auto* serializer_ptr = mock_serializer.get();

	env.solver = std::move(mock_solver);
	env.serializer = std::move(mock_serializer);

	env.run_single_timestep();
	EXPECT_TRUE(solver_ptr->solve_called);

	env.serialize_state(1.5);
	EXPECT_TRUE(serializer_ptr->serialize_called);
}

// Test set_mesh and has_mesh methods
TEST_F(EnvironmentTest, SetMeshAndHasMesh)
{
	environment env(0.1, 3, 2, 3);

	// Initially has_mesh returns false
	EXPECT_FALSE(env.has_mesh());

	// Set a mesh and verify has_mesh returns true
	auto mesh = physicore::cartesian_mesh(3, { -100, -100, -100 }, { 100, 100, 100 }, { 20, 20, 20 });
	env.set_mesh(mesh);
	EXPECT_TRUE(env.has_mesh());
}

// Test mesh operations and exception handling
TEST_F(EnvironmentTest, MeshOperationsAndErrorHandling)
{
	environment env(0.1, 3, 2, 3);

	// Test exception when getting mesh without setting it
	EXPECT_THROW(env.get_mesh(), std::runtime_error);

	// Set a mesh and verify retrieval
	auto mesh = physicore::cartesian_mesh(3, { -100, -100, -100 }, { 100, 100, 100 }, { 20, 20, 20 });
	env.set_mesh(mesh);

	const auto& retrieved = env.get_mesh();
	EXPECT_EQ(retrieved.dims, 3);
}

} // namespace physicore::mechanics::physicell::tests
