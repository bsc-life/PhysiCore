#include <filesystem>

#include <gtest/gtest.h>

#include "config_reader.h"

using namespace physicore::mechanics::physicell;

namespace {

std::filesystem::path settings_path()
{
	auto base = std::filesystem::path(__FILE__).parent_path() / "fixtures" / "PhysiCell_settings.xml";
	return std::filesystem::absolute(base);
}

} // namespace

TEST(ConfigReaderMechanics, ParsesMechanicsParametersFromSettings)
{
	auto config = parse_simulation_parameters(settings_path());
	const auto& params = config.parameters;

	ASSERT_EQ(config.cell_definition_names.size(), 3u);
	EXPECT_EQ(config.cell_definition_names[0], "director cell");
	ASSERT_EQ(config.substrate_names.size(), 2u);
	EXPECT_EQ(config.substrate_names[0], "director signal");

	// Domain use_2D propagates into params
	EXPECT_TRUE(params.is_2D);

	// Mechanics for director cell (ID 0)
	EXPECT_DOUBLE_EQ(params.cell_cell_adhesion_strength[0], 0.4);
	EXPECT_DOUBLE_EQ(params.cell_cell_repulsion_strength[0], 10.0);
	EXPECT_DOUBLE_EQ(params.relative_maximum_adhesion_distance[0], 1.25);
	EXPECT_DOUBLE_EQ(params.cell_adhesion_affinity[0][1], 1.0); // cargo cell affinity

	EXPECT_DOUBLE_EQ(params.attachment_elastic_coefficient[0], 0.5);
	EXPECT_DOUBLE_EQ(params.attachment_rate[0], 10.0);
	EXPECT_DOUBLE_EQ(params.detachment_rate[0], 0.0);

	// Motility for director cell
	EXPECT_DOUBLE_EQ(params.motility_speed[0], 1.0);
	EXPECT_DOUBLE_EQ(params.motility_persistence_time[0], 1.0);
	EXPECT_DOUBLE_EQ(params.motility_bias[0], 0.5);
	EXPECT_FALSE(params.is_movable[0]); // motility enabled is false in XML

	// Motility for cargo cell
	EXPECT_DOUBLE_EQ(params.motility_speed[1], 1.0);
	EXPECT_DOUBLE_EQ(params.motility_persistence_time[1], 1.0);
	EXPECT_DOUBLE_EQ(params.motility_bias[1], 0.5);
	EXPECT_FALSE(params.is_movable[1]); // motility enabled is false in XML


	// Chemotaxis matrices sized correctly and default disabled
	EXPECT_EQ(params.chemotaxis_sensitivity.size(), 3u);	// cell types
	EXPECT_EQ(params.chemotaxis_sensitivity[0].size(), 2u); // substrates
	EXPECT_FALSE(params.chemotaxis_enabled[0][0]);
	EXPECT_DOUBLE_EQ(params.chemotaxis_sensitivity[0][0], 0.0);
}
