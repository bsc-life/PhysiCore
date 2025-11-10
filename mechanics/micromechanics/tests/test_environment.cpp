#include <gtest/gtest.h>

#include "environment.h"

using namespace physicore::mechanics::micromechanics;

TEST(EnvironmentTest, RunSingleTimestep)
{
	environment env(0.1);
	env.run_single_timestep();
	// No assertion, just check it runs without error
}
