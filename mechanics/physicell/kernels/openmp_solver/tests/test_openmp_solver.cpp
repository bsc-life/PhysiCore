#include <gtest/gtest.h>

#include "physicell/environment.h"
#include "physicell/openmp_solver/register_solver.h"
#include "physicell/solver_registry.h"

using namespace physicore::mechanics::physicell;

TEST(OpenMPSolverTest, RegistersInRegistry)
{
	kernels::openmp_solver::attach_to_registry();

	auto solver = solver_registry::instance().get("openmp_solver");
	ASSERT_NE(solver, nullptr);
}

TEST(OpenMPSolverTest, CanRunViaEnvironment)
{
	kernels::openmp_solver::attach_to_registry();

	environment env(0.1);
	env.solver = solver_registry::instance().get("openmp_solver");

	ASSERT_NE(env.solver, nullptr);
	env.run_single_timestep();
}
