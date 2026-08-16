#include <gtest/gtest.h>

#include "physicell/solver.h"
#include "physicell/solver_registry.h"

namespace physicore::mechanics::physicell::tests {

using namespace physicore;
using namespace physicore::mechanics::physicell;

class RegistryTest : public ::testing::Test
{};

// Test registry covers the if body in registry.h::37 (key not found case)
TEST_F(RegistryTest, GetNonExistentSolverReturnsNullptr)
{
	// Try to get a solver that doesn't exist in the registry
	// This exercises the if (it == factories_.end()) branch at registry.h::37
	auto result = solver_registry::instance().get("nonexistent_solver_xyz");
	EXPECT_EQ(result, nullptr);
}

// Test registry covers the else case (key found)
TEST_F(RegistryTest, GetExistentSolverReturnsValidInstance)
{
	// Get the registered openmp solver
	// This exercises the else branch where the factory is found and called
	auto result = solver_registry::instance().get("openmp_solver");
	ASSERT_NE(result, nullptr);
	EXPECT_TRUE(dynamic_cast<solver*>(result.get()) != nullptr);
}

} // namespace physicore::mechanics::physicell::tests
