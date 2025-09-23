#include <gtest/gtest.h>

#include "solver_registry.h"

using namespace physicore;
using namespace physicore::biofvm;

TEST(SolverRegistryTest, CheckPresentSolvers)
{
	auto registry = solver_registry::instance();

	EXPECT_NE(registry.get("cpu_solver"), nullptr);
}
