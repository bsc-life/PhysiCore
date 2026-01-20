#include <common/base_agent_data.h>
#include <gtest/gtest.h>

#include "micromechanics/agent_container.h"
#include "micromechanics/agent_data.h"
#include "micromechanics/environment.h"
#include "micromechanics/uniform_grid_spatial_index.h"

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

TEST(UniformGridSpatialIndex, BuildAndQuery)
{
	// Setup environment
	environment env(0.01);

	auto base_data = std::make_unique<base_agent_data>();
	auto mech_data = std::make_unique<agent_data>(*base_data);

	env.agents = std::make_unique<agent_container>(std::move(base_data), std::move(mech_data));

	// Add some agents
	auto& agents = *env.agents;
	auto& mech_data_ref = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);
	auto& base_data_ref = *std::get<std::unique_ptr<base_agent_data>>(agents.agent_datas);

	// Agent 0 at (0,0,0)
	agents.create();
	base_data_ref.positions[0] = 0.0;
	base_data_ref.positions[1] = 0.0;
	base_data_ref.positions[2] = 0.0;
	mech_data_ref.radii[0] = 10.0;

	// Agent 1 at (15,0,0) - should be neighbor of 0 (dist 15 < 20)
	agents.create();
	base_data_ref.positions[3] = 15.0;
	base_data_ref.positions[4] = 0.0;
	base_data_ref.positions[5] = 0.0;
	mech_data_ref.radii[1] = 10.0;

	// Agent 2 at (100,0,0) - far away
	agents.create();
	base_data_ref.positions[6] = 100.0;
	base_data_ref.positions[7] = 0.0;
	base_data_ref.positions[8] = 0.0;
	mech_data_ref.radii[2] = 10.0;

	// Build index
	uniform_grid_spatial_index index;
	index.build(env);

	// Query neighbors for Agent 0 with radius 20
	auto neighbors = index.query_neighbors(env, 0, 20.0);

	bool found_1 = false;
	bool found_2 = false;

	for (auto idx : neighbors)
	{
		if (idx == 1)
			found_1 = true;
		if (idx == 2)
			found_2 = true;
	}

	EXPECT_TRUE(found_1);
	EXPECT_FALSE(found_2);
}
