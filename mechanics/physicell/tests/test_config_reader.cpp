#include <filesystem>
#include <fstream>
#include <stdexcept>

#include <common/types.h>
#include <gtest/gtest.h>

#include "config_reader.h"

using namespace physicore;
using namespace physicore::mechanics::physicell;

// Helper to create a complete test XML with domain and overall sections
static std::string create_complete_xml(const std::string& microenv_section, const std::string& cell_defs_section,
									   const std::string& domain_section = "", const std::string& overall_section = "")
{
	const std::string domain = domain_section.empty() ? R"(
<domain>
<x_min>-500</x_min>
<x_max>500</x_max>
<y_min>-500</y_min>
<y_max>500</y_max>
<z_min>-10</z_min>
<z_max>10</z_max>
<dx>20</dx>
<dy>20</dy>
<dz>20</dz>
<use_2D>true</use_2D>
</domain>)"
													  : domain_section;

	const std::string overall = overall_section.empty() ? R"(
<overall>
<max_time>14400</max_time>
<time_units>min</time_units>
<space_units>micron</space_units>
<dt_mechanics>0.1</dt_mechanics>
</overall>)"
														: overall_section;

	return R"(<?xml version="1.0"?>
<PhysiCell_settings version="devel-version">)"
		   + domain + overall + microenv_section + cell_defs_section + R"(
</PhysiCell_settings>)";
}

class MechanicsConfigReaderTest : public ::testing::Test
{
protected:
	std::filesystem::path test_config_file;

	void SetUp() override { test_config_file = "test_mechanics_config.xml"; }

	void TearDown() override
	{
		if (std::filesystem::exists(test_config_file))
		{
			std::filesystem::remove(test_config_file);
		}
	}

	void write_config(const std::string& xml_content)
	{
		std::ofstream ofs(test_config_file);
		ofs << xml_content;
		ofs.close();
	}
};

// ============================================================================
// TESTS FOR DOMAIN CONFIGURATION
// ============================================================================

TEST_F(MechanicsConfigReaderTest, DomainConfig_ParsesAllFields)
{
	const std::string custom_domain = R"(
<domain>
<x_min>-1000</x_min>
<x_max>1000</x_max>
<y_min>-750</y_min>
<y_max>750</y_max>
<z_min>-50</z_min>
<z_max>50</z_max>
<dx>25</dx>
<dy>25</dy>
<dz>25</dz>
<use_2D>false</use_2D>
</domain>)";

	const std::string microenv = R"(
<microenvironment_setup>
<variable name="oxygen" ID="0" />
</microenvironment_setup>)";

	const std::string cells = R"(
<cell_definitions>
<cell_definition name="default cell" ID="0">
<phenotype>
<mechanics>
<cell_cell_adhesion_strength>0.4</cell_cell_adhesion_strength>
<cell_cell_repulsion_strength>10.0</cell_cell_repulsion_strength>
<relative_maximum_adhesion_distance>1.25</relative_maximum_adhesion_distance>
<cell_adhesion_affinities>
<cell_adhesion_affinity name="default cell">1.0</cell_adhesion_affinity>
</cell_adhesion_affinities>
<attachment_elastic_constant>0.5</attachment_elastic_constant>
<attachment_rate>10.0</attachment_rate>
<detachment_rate>0.1</detachment_rate>
</mechanics>
<motility>
<speed>2.0</speed>
<persistence_time>5.0</persistence_time>
<migration_bias>0.5</migration_bias>
<options>
<enabled>false</enabled>
<use_2D>false</use_2D>
</options>
</motility>
</phenotype>
</cell_definition>
</cell_definitions>)";

	write_config(create_complete_xml(microenv, cells, custom_domain));
	auto config = parse_simulation_parameters(test_config_file);

	// Verify domain configuration
	EXPECT_DOUBLE_EQ(config.domain.x_min, -1000.0);
	EXPECT_DOUBLE_EQ(config.domain.x_max, 1000.0);
	EXPECT_DOUBLE_EQ(config.domain.y_min, -750.0);
	EXPECT_DOUBLE_EQ(config.domain.y_max, 750.0);
	EXPECT_DOUBLE_EQ(config.domain.z_min, -50.0);
	EXPECT_DOUBLE_EQ(config.domain.z_max, 50.0);
	EXPECT_DOUBLE_EQ(config.domain.dx, 25.0);
	EXPECT_DOUBLE_EQ(config.domain.dy, 25.0);
	EXPECT_DOUBLE_EQ(config.domain.dz, 25.0);
	EXPECT_FALSE(config.domain.use_2D);
	EXPECT_FALSE(config.is_2D); // Should match domain.use_2D
}

TEST_F(MechanicsConfigReaderTest, DomainConfig_2DFlag)
{
	const std::string custom_domain = R"(
<domain>
<x_min>-500</x_min>
<x_max>500</x_max>
<y_min>-500</y_min>
<y_max>500</y_max>
<z_min>-10</z_min>
<z_max>10</z_max>
<dx>20</dx>
<dy>20</dy>
<dz>20</dz>
<use_2D>true</use_2D>
</domain>)";

	const std::string microenv = R"(
<microenvironment_setup>
<variable name="oxygen" ID="0" />
</microenvironment_setup>)";

	const std::string cells = R"(
<cell_definitions>
<cell_definition name="cell1" ID="0">
<phenotype>
<mechanics>
<cell_cell_adhesion_strength>0.0</cell_cell_adhesion_strength>
<cell_cell_repulsion_strength>0.0</cell_cell_repulsion_strength>
<relative_maximum_adhesion_distance>1.0</relative_maximum_adhesion_distance>
<cell_adhesion_affinities>
<cell_adhesion_affinity name="cell1">1.0</cell_adhesion_affinity>
</cell_adhesion_affinities>
<attachment_elastic_constant>0.0</attachment_elastic_constant>
<attachment_rate>0.0</attachment_rate>
<detachment_rate>0.0</detachment_rate>
</mechanics>
<motility>
<speed>0</speed>
<persistence_time>0</persistence_time>
<migration_bias>0</migration_bias>
<options>
<enabled>false</enabled>
<use_2D>true</use_2D>
</options>
</motility>
</phenotype>
</cell_definition>
</cell_definitions>)";

	write_config(create_complete_xml(microenv, cells, custom_domain));
	auto config = parse_simulation_parameters(test_config_file);

	EXPECT_TRUE(config.domain.use_2D);
	EXPECT_TRUE(config.is_2D);
}

// ============================================================================
// TESTS FOR OVERALL CONFIGURATION
// ============================================================================

TEST_F(MechanicsConfigReaderTest, OverallConfig_ParsesAllFields)
{
	const std::string custom_overall = R"(
<overall>
<max_time>28800</max_time>
<time_units>hours</time_units>
<space_units>millimeters</space_units>
<dt_mechanics>0.05</dt_mechanics>
</overall>)";

	const std::string microenv = R"(
<microenvironment_setup>
<variable name="glucose" ID="0" />
</microenvironment_setup>)";

	const std::string cells = R"(
<cell_definitions>
<cell_definition name="default" ID="0">
<phenotype>
<mechanics>
<cell_cell_adhesion_strength>0.0</cell_cell_adhesion_strength>
<cell_cell_repulsion_strength>0.0</cell_cell_repulsion_strength>
<relative_maximum_adhesion_distance>1.0</relative_maximum_adhesion_distance>
<cell_adhesion_affinities>
<cell_adhesion_affinity name="default">1.0</cell_adhesion_affinity>
</cell_adhesion_affinities>
<attachment_elastic_constant>0.0</attachment_elastic_constant>
<attachment_rate>0.0</attachment_rate>
<detachment_rate>0.0</detachment_rate>
</mechanics>
<motility>
<speed>0</speed>
<persistence_time>0</persistence_time>
<migration_bias>0</migration_bias>
<options>
<enabled>false</enabled>
<use_2D>false</use_2D>
</options>
</motility>
</phenotype>
</cell_definition>
</cell_definitions>)";

	write_config(create_complete_xml(microenv, cells, "", custom_overall));
	auto config = parse_simulation_parameters(test_config_file);

	// Verify overall configuration
	EXPECT_DOUBLE_EQ(config.overall.max_time, 28800.0);
	EXPECT_EQ(config.overall.time_units, "hours");
	EXPECT_EQ(config.overall.space_units, "millimeters");
	EXPECT_DOUBLE_EQ(config.overall.dt_mechanics, 0.05);
}

TEST_F(MechanicsConfigReaderTest, OverallConfig_DifferentTimeUnits)
{
	const std::string custom_overall = R"(
<overall>
<max_time>3600</max_time>
<time_units>seconds</time_units>
<space_units>micron</space_units>
<dt_mechanics>0.1</dt_mechanics>
</overall>)";

	const std::string microenv = R"(
<microenvironment_setup>
<variable name="substrate1" ID="0" />
</microenvironment_setup>)";

	const std::string cells = R"(
<cell_definitions>
<cell_definition name="cell" ID="0">
<phenotype>
<mechanics>
<cell_cell_adhesion_strength>0.0</cell_cell_adhesion_strength>
<cell_cell_repulsion_strength>0.0</cell_cell_repulsion_strength>
<relative_maximum_adhesion_distance>1.0</relative_maximum_adhesion_distance>
<cell_adhesion_affinities>
<cell_adhesion_affinity name="cell">1.0</cell_adhesion_affinity>
</cell_adhesion_affinities>
<attachment_elastic_constant>0.0</attachment_elastic_constant>
<attachment_rate>0.0</attachment_rate>
<detachment_rate>0.0</detachment_rate>
</mechanics>
<motility>
<speed>0</speed>
<persistence_time>0</persistence_time>
<migration_bias>0</migration_bias>
<options>
<enabled>false</enabled>
<use_2D>false</use_2D>
</options>
</motility>
</phenotype>
</cell_definition>
</cell_definitions>)";

	write_config(create_complete_xml(microenv, cells, "", custom_overall));
	auto config = parse_simulation_parameters(test_config_file);

	EXPECT_EQ(config.overall.time_units, "seconds");
	EXPECT_EQ(config.overall.space_units, "micron");
}

// ============================================================================
// POSITIVE TESTS FOR MECHANICS PARAMETERS
// ============================================================================

TEST_F(MechanicsConfigReaderTest, ParsesBasicMechanicsParameters)
{
	const std::string microenv = R"(
<microenvironment_setup>
<variable name="oxygen" ID="0" />
<variable name="glucose" ID="1" />
</microenvironment_setup>)";

	const std::string cells = R"(
<cell_definitions>
<cell_definition name="default cell" ID="0">
<phenotype>
<mechanics>
<cell_cell_adhesion_strength>0.4</cell_cell_adhesion_strength>
<cell_cell_repulsion_strength>10.0</cell_cell_repulsion_strength>
<relative_maximum_adhesion_distance>1.25</relative_maximum_adhesion_distance>
<cell_adhesion_affinities>
<cell_adhesion_affinity name="default cell">1.0</cell_adhesion_affinity>
</cell_adhesion_affinities>
<attachment_elastic_constant>0.5</attachment_elastic_constant>
<attachment_rate>10.0</attachment_rate>
<detachment_rate>0.1</detachment_rate>
</mechanics>
<motility>
<speed>2.0</speed>
<persistence_time>5.0</persistence_time>
<migration_bias>0.5</migration_bias>
<options>
<enabled>true</enabled>
<use_2D>true</use_2D>
</options>
</motility>
</phenotype>
</cell_definition>
</cell_definitions>)";

	write_config(create_complete_xml(microenv, cells));

	auto config = parse_simulation_parameters(test_config_file);

	// Basic verification
	ASSERT_EQ(config.cell_types.size(), 1u);
	EXPECT_EQ(config.cell_types[0].name, "default cell");

	// Mechanics parameters
	const auto& params = config.cell_types[0];
	EXPECT_DOUBLE_EQ(params.cell_cell_adhesion_strength, 0.4);
	EXPECT_DOUBLE_EQ(params.cell_cell_repulsion_strength, 10.0);
	EXPECT_DOUBLE_EQ(params.relative_maximum_adhesion_distance, 1.25);

	// Attachment parameters
	EXPECT_DOUBLE_EQ(params.attachment_elastic_coefficient, 0.5);
	EXPECT_DOUBLE_EQ(params.attachment_rate, 10.0);
	EXPECT_DOUBLE_EQ(params.detachment_rate, 0.1);

	// Motility parameters
	EXPECT_TRUE(params.is_movable);
	EXPECT_DOUBLE_EQ(params.motility_speed, 2.0);
	EXPECT_DOUBLE_EQ(params.motility_persistence_time, 5.0);
	EXPECT_DOUBLE_EQ(params.motility_bias, 0.5);
}

TEST_F(MechanicsConfigReaderTest, CornerCase_DefaultDomain2DFalse)
{
	// When domain section is missing, should default to use_2D = false
	const std::string no_domain = R"(
<overall>
<max_time>14400</max_time>
<time_units>min</time_units>
<space_units>micron</space_units>
<dt_mechanics>0.1</dt_mechanics>
</overall>)";

	const std::string microenv = R"(
<microenvironment_setup>
<variable name="v1" ID="0" />
</microenvironment_setup>)";

	const std::string cells = R"(
<cell_definitions>
<cell_definition name="default" ID="0">
<phenotype>
<mechanics>
<cell_cell_adhesion_strength>0.0</cell_cell_adhesion_strength>
<cell_cell_repulsion_strength>0.0</cell_cell_repulsion_strength>
<relative_maximum_adhesion_distance>1.0</relative_maximum_adhesion_distance>
<cell_adhesion_affinities>
<cell_adhesion_affinity name="default">1.0</cell_adhesion_affinity>
</cell_adhesion_affinities>
<attachment_elastic_constant>0.0</attachment_elastic_constant>
<attachment_rate>0.0</attachment_rate>
<detachment_rate>0.0</detachment_rate>
</mechanics>
<motility>
<speed>0</speed>
<persistence_time>0</persistence_time>
<migration_bias>0</migration_bias>
<options>
<enabled>false</enabled>
<use_2D>false</use_2D>
</options>
</motility>
</phenotype>
</cell_definition>
</cell_definitions>)";

	const std::string xml = R"(<?xml version="1.0"?>
<PhysiCell_settings version="devel-version">)"
							+ no_domain + microenv + cells + R"(
</PhysiCell_settings>)";

	write_config(xml);
	auto config = parse_simulation_parameters(test_config_file);

	EXPECT_FALSE(config.is_2D);
}

TEST_F(MechanicsConfigReaderTest, CornerCase_NoSubstrates)
{
	const std::string microenv = R"(
<microenvironment_setup>
</microenvironment_setup>)";

	const std::string cells = R"(
<cell_definitions>
<cell_definition name="cell1" ID="0">
<phenotype>
<mechanics>
<cell_cell_adhesion_strength>0.0</cell_cell_adhesion_strength>
<cell_cell_repulsion_strength>0.0</cell_cell_repulsion_strength>
<relative_maximum_adhesion_distance>1.0</relative_maximum_adhesion_distance>
<cell_adhesion_affinities>
<cell_adhesion_affinity name="cell1">1.0</cell_adhesion_affinity>
</cell_adhesion_affinities>
<attachment_elastic_constant>0.0</attachment_elastic_constant>
<attachment_rate>0.0</attachment_rate>
<detachment_rate>0.0</detachment_rate>
</mechanics>
<motility>
<speed>0</speed>
<persistence_time>0</persistence_time>
<migration_bias>0</migration_bias>
<options>
<enabled>false</enabled>
<use_2D>true</use_2D>
</options>
</motility>
</phenotype>
</cell_definition>
</cell_definitions>)";

	write_config(create_complete_xml(microenv, cells));
	auto config = parse_simulation_parameters(test_config_file);

	// Config should parse successfully even without substrates
	ASSERT_EQ(config.cell_types.size(), 1u);
}

// ============================================================================
// NEGATIVE TESTS
// ============================================================================

TEST_F(MechanicsConfigReaderTest, Negative_MissingRequiredElement)
{
	const std::string microenv = R"(
<microenvironment_setup>
<variable name="oxygen" ID="0" />
</microenvironment_setup>)";

	const std::string cells = R"(
<cell_definitions>
<cell_definition name="cell1" ID="0">
<phenotype>
<mechanics>
<!-- Missing cell_cell_adhesion_strength (required) -->
<cell_cell_repulsion_strength>0.0</cell_cell_repulsion_strength>
<relative_maximum_adhesion_distance>1.0</relative_maximum_adhesion_distance>
<cell_adhesion_affinities>
<cell_adhesion_affinity name="cell1">1.0</cell_adhesion_affinity>
</cell_adhesion_affinities>
<attachment_elastic_constant>0.0</attachment_elastic_constant>
<attachment_rate>0.0</attachment_rate>
<detachment_rate>0.0</detachment_rate>
</mechanics>
<motility>
<speed>0</speed>
<persistence_time>0</persistence_time>
<migration_bias>0</migration_bias>
<options>
<enabled>false</enabled>
<use_2D>true</use_2D>
</options>
</motility>
</phenotype>
</cell_definition>
</cell_definitions>)";

	write_config(create_complete_xml(microenv, cells));

	// Missing required element should throw
	EXPECT_THROW(parse_simulation_parameters(test_config_file), std::runtime_error);
}

TEST_F(MechanicsConfigReaderTest, Negative_UnknownCellTypeInAffinity)
{
	const std::string microenv = R"(
<microenvironment_setup>
<variable name="oxygen" ID="0" />
</microenvironment_setup>)";

	const std::string cells = R"(
<cell_definitions>
<cell_definition name="cell1" ID="0">
<phenotype>
<mechanics>
<cell_cell_adhesion_strength>0.0</cell_cell_adhesion_strength>
<cell_cell_repulsion_strength>0.0</cell_cell_repulsion_strength>
<relative_maximum_adhesion_distance>1.0</relative_maximum_adhesion_distance>
<cell_adhesion_affinities>
<cell_adhesion_affinity name="cell1">1.0</cell_adhesion_affinity>
<cell_adhesion_affinity name="unknown_cell">0.5</cell_adhesion_affinity>
</cell_adhesion_affinities>
<attachment_elastic_constant>0.0</attachment_elastic_constant>
<attachment_rate>0.0</attachment_rate>
<detachment_rate>0.0</detachment_rate>
</mechanics>
<motility>
<speed>0</speed>
<persistence_time>0</persistence_time>
<migration_bias>0</migration_bias>
<options>
<enabled>false</enabled>
<use_2D>true</use_2D>
</options>
</motility>
</phenotype>
</cell_definition>
</cell_definitions>)";

	write_config(create_complete_xml(microenv, cells));

	// Unknown cell type should throw
	EXPECT_THROW(parse_simulation_parameters(test_config_file), std::runtime_error);
}
