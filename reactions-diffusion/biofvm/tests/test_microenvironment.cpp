#include <gtest/gtest.h>

#include "microenvironment.h"

using namespace physicore::biofvm;

TEST(MicroenvironmentTest, RunSingleTimestep)
{
	microenvironment env(0.1);
	env.run_single_timestep();
	// No assertion, just check it runs without error
}
