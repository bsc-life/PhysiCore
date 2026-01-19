#include <memory>

#include <common/base_agent_data.h>
#include <gtest/gtest.h>

#include "micromechanics/agent_container.h"
#include "micromechanics/agent_data.h"
#include "micromechanics/environment.h"
#include "micromechanics/solver_registry.h"
#include "micromechanics/uniform_grid_spatial_index.h"

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

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
