#include <string>

#include <common/types.h>
#include <gtest/gtest.h>

#include "physicell/solver_registry.h"

using namespace physicore;
using namespace physicore::mechanics::physicell;

TEST(SolverRegistryTest, OpenMPSolverIsRegistered)
{
	auto solver_instance = solver_registry::instance().get("openmp_solver");
	ASSERT_NE(solver_instance, nullptr);
}
