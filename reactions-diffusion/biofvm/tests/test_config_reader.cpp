#include <filesystem>
#include <fstream>

#include <gtest/gtest.h>

#include "config_reader.h"
#include "microenvironment.h"

using namespace physicore;
using namespace physicore::biofvm;

class ConfigReaderTest : public ::testing::Test
{
protected:
	std::filesystem::path test_config_file;

	void SetUp() override
	{
		// Create a test config file with known values
		test_config_file = "test_config.xml";
		std::ofstream ofs(test_config_file);
		ofs << R"(<?xml version="1.0"?>
<PhysiCell_settings version="devel-version">
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
	</domain>

	<overall>
		<max_time units="min">14400</max_time>
		<time_units>min</time_units>
		<space_units>micron</space_units>
		<dt_diffusion units="min">0.01</dt_diffusion>
		<dt_mechanics units="min">0.1</dt_mechanics>
		<dt_phenotype units="min">6</dt_phenotype>
	</overall>

	<microenvironment_setup>
		<variable name="oxygen" units="dimensionless" ID="0">
			<physical_parameter_set>
				<diffusion_coefficient units="micron^2/min">100000.0</diffusion_coefficient>
				<decay_rate units="1/min">0.1</decay_rate>
			</physical_parameter_set>
			<initial_condition units="mmHg">38</initial_condition>
			<Dirichlet_options>
				<boundary_value ID="xmin" enabled="True">38</boundary_value>
				<boundary_value ID="xmax" enabled="True">10</boundary_value>
				<boundary_value ID="ymin" enabled="True">10</boundary_value>
				<boundary_value ID="ymax" enabled="True">38</boundary_value>
				<boundary_value ID="zmin" enabled="False">0</boundary_value>
				<boundary_value ID="zmax" enabled="False">0</boundary_value>
			</Dirichlet_options>
		</variable>
		<variable name="necrotic debris" units="dimensionless" ID="1">
			<physical_parameter_set>
				<diffusion_coefficient units="micron^2/min">10</diffusion_coefficient>
				<decay_rate units="1/min">0.1</decay_rate>
			</physical_parameter_set>
			<initial_condition units="mmHg">0</initial_condition>
			<Dirichlet_options>
				<boundary_value ID="xmin" enabled="False">0</boundary_value>
				<boundary_value ID="xmax" enabled="False">0</boundary_value>
				<boundary_value ID="ymin" enabled="False">0</boundary_value>
				<boundary_value ID="ymax" enabled="False">0</boundary_value>
				<boundary_value ID="zmin" enabled="False">0</boundary_value>
				<boundary_value ID="zmax" enabled="False">0</boundary_value>
			</Dirichlet_options>
		</variable>
		<variable name="apoptotic debris" units="dimensionless" ID="2">
			<physical_parameter_set>
				<diffusion_coefficient units="micron^2/min">10</diffusion_coefficient>
				<decay_rate units="1/min">0.1</decay_rate>
			</physical_parameter_set>
			<initial_condition units="mmHg">0</initial_condition>
			<Dirichlet_options>
				<boundary_value ID="xmin" enabled="False">0</boundary_value>
				<boundary_value ID="xmax" enabled="False">0</boundary_value>
				<boundary_value ID="ymin" enabled="False">0</boundary_value>
				<boundary_value ID="ymax" enabled="False">0</boundary_value>
				<boundary_value ID="zmin" enabled="False">0</boundary_value>
				<boundary_value ID="zmax" enabled="False">0</boundary_value>
			</Dirichlet_options>
		</variable>
		<options>
			<calculate_gradients>true</calculate_gradients>
			<track_internalized_substrates_in_each_agent>true</track_internalized_substrates_in_each_agent>
		</options>
	</microenvironment_setup>
</PhysiCell_settings>
)";
		ofs.close();
	}

	void TearDown() override
	{
		// Clean up test config file
		if (std::filesystem::exists(test_config_file))
		{
			std::filesystem::remove(test_config_file);
		}
	}
};

TEST_F(ConfigReaderTest, ParsePhysiCellConfigFile)
{
	physicell_config config = parse_physicell_config(test_config_file);

	// Verify domain configuration
	EXPECT_EQ(config.domain.x_min, -500.0);
	EXPECT_EQ(config.domain.x_max, 500.0);
	EXPECT_EQ(config.domain.y_min, -500.0);
	EXPECT_EQ(config.domain.y_max, 500.0);
	EXPECT_EQ(config.domain.z_min, -10.0);
	EXPECT_EQ(config.domain.z_max, 10.0);
	EXPECT_EQ(config.domain.dx, 20.0);
	EXPECT_EQ(config.domain.dy, 20.0);
	EXPECT_EQ(config.domain.dz, 20.0);
	EXPECT_TRUE(config.domain.use_2D);

	// Verify overall configuration
	EXPECT_EQ(config.overall.max_time, 14400.0);
	EXPECT_EQ(config.overall.time_units, "min");
	EXPECT_EQ(config.overall.space_units, "micron");
	EXPECT_EQ(config.overall.dt_diffusion, 0.01);
	EXPECT_EQ(config.overall.dt_mechanics, 0.1);
	EXPECT_EQ(config.overall.dt_phenotype, 6.0);

	// Verify microenvironment configuration
	EXPECT_TRUE(config.microenvironment.calculate_gradients);
	EXPECT_TRUE(config.microenvironment.track_internalized_substrates);

	// Should have 3 substrates: oxygen, necrotic debris, apoptotic debris
	ASSERT_GE(config.microenvironment.variables.size(), 3);

	// Verify oxygen substrate (first variable)
	const auto& oxygen = config.microenvironment.variables[0];
	EXPECT_EQ(oxygen.name, "oxygen");
	EXPECT_EQ(oxygen.units, "dimensionless");
	EXPECT_EQ(oxygen.id, 0);
	EXPECT_EQ(oxygen.diffusion_coefficient, 100000.0);
	EXPECT_EQ(oxygen.decay_rate, 0.1);
	EXPECT_EQ(oxygen.initial_condition, 38.0);

	// Verify oxygen Dirichlet boundary conditions
	EXPECT_EQ(oxygen.boundary_conditions.mins_values[0], 38.0); // xmin
	EXPECT_EQ(oxygen.boundary_conditions.maxs_values[0], 10.0); // xmax
	EXPECT_EQ(oxygen.boundary_conditions.mins_values[1], 10.0); // ymin
	EXPECT_EQ(oxygen.boundary_conditions.maxs_values[1], 38.0); // ymax
	EXPECT_EQ(oxygen.boundary_conditions.mins_values[2], 0.0);	// zmin
	EXPECT_EQ(oxygen.boundary_conditions.maxs_values[2], 0.0);	// zmax

	EXPECT_TRUE(oxygen.boundary_conditions.mins_conditions[0]);	 // xmin enabled
	EXPECT_TRUE(oxygen.boundary_conditions.maxs_conditions[0]);	 // xmax enabled
	EXPECT_TRUE(oxygen.boundary_conditions.mins_conditions[1]);	 // ymin enabled
	EXPECT_TRUE(oxygen.boundary_conditions.maxs_conditions[1]);	 // ymax enabled
	EXPECT_FALSE(oxygen.boundary_conditions.mins_conditions[2]); // zmin disabled
	EXPECT_FALSE(oxygen.boundary_conditions.maxs_conditions[2]); // zmax disabled
}

TEST_F(ConfigReaderTest, MissingFileThrowsException)
{
	const std::filesystem::path nonexistent = "nonexistent_config.xml";
	EXPECT_THROW(parse_physicell_config(nonexistent), std::runtime_error);
}

TEST_F(ConfigReaderTest, MalformedXMLThrowsException)
{
	// Create a temporary malformed XML file
	const std::filesystem::path temp_file = "malformed_test.xml";
	std::ofstream ofs(temp_file);
	ofs << "<PhysiCell_settings>\n<domain>\n<x_min>100</x_min>\n"; // Missing closing tags
	ofs.close();

	EXPECT_THROW(parse_physicell_config(temp_file), std::runtime_error);

	// Cleanup
	std::filesystem::remove(temp_file);
}

TEST_F(ConfigReaderTest, MissingRequiredTagThrowsException)
{
	// Create a temporary XML file missing required <overall> tag
	const std::filesystem::path temp_file = "missing_tag_test.xml";
	std::ofstream ofs(temp_file);
	ofs << R"(<?xml version="1.0"?>
<PhysiCell_settings>
	<domain>
		<x_min>0</x_min>
		<x_max>100</x_max>
		<y_min>0</y_min>
		<y_max>100</y_max>
		<z_min>0</z_min>
		<z_max>100</z_max>
		<dx>10</dx>
		<dy>10</dy>
		<dz>10</dz>
		<use_2D>false</use_2D>
	</domain>
	<!-- Missing <overall> tag -->
	<microenvironment_setup>
		<variable name="test" units="dimensionless" ID="0">
			<physical_parameter_set>
				<diffusion_coefficient units="micron^2/min">100.0</diffusion_coefficient>
				<decay_rate units="1/min">0.1</decay_rate>
			</physical_parameter_set>
			<initial_condition units="mmHg">1.0</initial_condition>
		</variable>
		<options>
			<calculate_gradients>false</calculate_gradients>
			<track_internalized_substrates_in_each_agent>false</track_internalized_substrates_in_each_agent>
		</options>
	</microenvironment_setup>
</PhysiCell_settings>
)";
	ofs.close();

	EXPECT_THROW(parse_physicell_config(temp_file), std::runtime_error);

	// Cleanup
	std::filesystem::remove(temp_file);
}

TEST_F(ConfigReaderTest, CreateMicroenvironmentFromConfig)
{
	// Create microenvironment from config file
	auto microenv = microenvironment::create_from_config(test_config_file);

	ASSERT_NE(microenv, nullptr);

	// Verify basic properties
	EXPECT_EQ(microenv->time_units, "min");
	EXPECT_EQ(microenv->space_units, "micron");
	EXPECT_EQ(microenv->diffusion_timestep, 0.01);

	// Verify mesh configuration
	EXPECT_EQ(microenv->mesh.dims, 2); // use_2D is true
	EXPECT_EQ(microenv->mesh.bounding_box_mins[0], -500);
	EXPECT_EQ(microenv->mesh.bounding_box_maxs[0], 500);
	EXPECT_EQ(microenv->mesh.voxel_shape[0], 20);
	EXPECT_EQ(microenv->mesh.voxel_shape[1], 20);

	// Verify substrates
	EXPECT_GE(microenv->substrates_count, 3);
	EXPECT_EQ(microenv->substrates_names[0], "oxygen");

	// Verify diffusion coefficient and decay rate for oxygen
	EXPECT_EQ(microenv->diffusion_coefficients[0], 100000.0);
	EXPECT_EQ(microenv->decay_rates[0], 0.1);
	EXPECT_EQ(microenv->initial_conditions[0], 38.0);

	// Verify compute internalized substrates option
	EXPECT_TRUE(microenv->compute_internalized_substrates);

	// Verify Dirichlet boundary conditions exist
	ASSERT_NE(microenv->dirichlet_min_boundary_values[0], nullptr);
	ASSERT_NE(microenv->dirichlet_max_boundary_values[0], nullptr);
	ASSERT_NE(microenv->dirichlet_min_boundary_conditions[0], nullptr);
	ASSERT_NE(microenv->dirichlet_max_boundary_conditions[0], nullptr);

	// Verify oxygen boundary values
	EXPECT_EQ(microenv->dirichlet_min_boundary_values[0][0], 38.0); // xmin
	EXPECT_EQ(microenv->dirichlet_max_boundary_values[0][0], 10.0); // xmax
	EXPECT_TRUE(microenv->dirichlet_min_boundary_conditions[0][0]); // xmin enabled
	EXPECT_TRUE(microenv->dirichlet_max_boundary_conditions[0][0]); // xmax enabled
}

TEST(ConfigReaderErrorTest, MissingRootNodeThrowsException)
{
	// Create XML without PhysiCell_settings root
	const std::filesystem::path temp_file = "missing_root_test.xml";
	std::ofstream ofs(temp_file);
	ofs << R"(<?xml version="1.0"?>
<Settings>
	<domain>
		<x_min>0</x_min>
		<x_max>100</x_max>
	</domain>
</Settings>
)";
	ofs.close();

	EXPECT_THROW(parse_physicell_config(temp_file), std::runtime_error);

	std::filesystem::remove(temp_file);
}

TEST(ConfigReaderErrorTest, MissingOptionsNodeUsesDefaults)
{
	// Create XML with microenvironment but no options
	const std::filesystem::path temp_file = "missing_options_test.xml";
	std::ofstream ofs(temp_file);
	ofs << R"(<?xml version="1.0"?>
<PhysiCell_settings>
	<domain>
		<x_min>0</x_min>
		<x_max>100</x_max>
		<y_min>0</y_min>
		<y_max>100</y_max>
		<z_min>0</z_min>
		<z_max>100</z_max>
		<dx>10</dx>
		<dy>10</dy>
		<dz>10</dz>
		<use_2D>false</use_2D>
	</domain>

	<overall>
		<max_time units="min">100</max_time>
		<time_units>min</time_units>
		<space_units>micron</space_units>
		<dt_diffusion units="min">0.01</dt_diffusion>
		<dt_mechanics units="min">0.1</dt_mechanics>
		<dt_phenotype units="min">6</dt_phenotype>
	</overall>

	<microenvironment_setup>
		<variable name="test" units="dimensionless" ID="0">
			<physical_parameter_set>
				<diffusion_coefficient units="micron^2/min">100.0</diffusion_coefficient>
				<decay_rate units="1/min">0.1</decay_rate>
			</physical_parameter_set>
			<initial_condition units="mmHg">1.0</initial_condition>
		</variable>
	</microenvironment_setup>
</PhysiCell_settings>
)";
	ofs.close();

	const physicell_config config = parse_physicell_config(temp_file);

	// Options should default to false
	EXPECT_FALSE(config.microenvironment.calculate_gradients);
	EXPECT_FALSE(config.microenvironment.track_internalized_substrates);

	std::filesystem::remove(temp_file);
}

TEST(ConfigReaderErrorTest, NoSubstratesThrowsException)
{
	// Create XML with microenvironment_setup but no variables
	const std::filesystem::path temp_file = "no_substrates_test.xml";
	std::ofstream ofs(temp_file);
	ofs << R"(<?xml version="1.0"?>
<PhysiCell_settings>
	<domain>
		<x_min>0</x_min>
		<x_max>100</x_max>
		<y_min>0</y_min>
		<y_max>100</y_max>
		<z_min>0</z_min>
		<z_max>100</z_max>
		<dx>10</dx>
		<dy>10</dy>
		<dz>10</dz>
		<use_2D>false</use_2D>
	</domain>

	<overall>
		<max_time units="min">100</max_time>
		<time_units>min</time_units>
		<space_units>micron</space_units>
		<dt_diffusion units="min">0.01</dt_diffusion>
		<dt_mechanics units="min">0.1</dt_mechanics>
		<dt_phenotype units="min">6</dt_phenotype>
	</overall>

	<microenvironment_setup>
		<options>
			<calculate_gradients>false</calculate_gradients>
			<track_internalized_substrates_in_each_agent>false</track_internalized_substrates_in_each_agent>
		</options>
	</microenvironment_setup>
</PhysiCell_settings>
)";
	ofs.close();

	EXPECT_THROW(parse_physicell_config(temp_file), std::runtime_error);

	std::filesystem::remove(temp_file);
}

TEST(ConfigReaderErrorTest, LegacyDirichletBoundaryCondition)
{
	// Create XML using legacy Dirichlet_boundary_condition instead of Dirichlet_options
	const std::filesystem::path temp_file = "legacy_dirichlet_test.xml";
	std::ofstream ofs(temp_file);
	ofs << R"(<?xml version="1.0"?>
<PhysiCell_settings>
	<domain>
		<x_min>-100</x_min>
		<x_max>100</x_max>
		<y_min>-100</y_min>
		<y_max>100</y_max>
		<z_min>-100</z_min>
		<z_max>100</z_max>
		<dx>10</dx>
		<dy>10</dy>
		<dz>10</dz>
		<use_2D>false</use_2D>
	</domain>

	<overall>
		<max_time units="min">100</max_time>
		<time_units>min</time_units>
		<space_units>micron</space_units>
		<dt_diffusion units="min">0.01</dt_diffusion>
		<dt_mechanics units="min">0.1</dt_mechanics>
		<dt_phenotype units="min">6</dt_phenotype>
	</overall>

	<microenvironment_setup>
		<variable name="oxygen" units="mmHg" ID="0">
			<physical_parameter_set>
				<diffusion_coefficient units="micron^2/min">100000.0</diffusion_coefficient>
				<decay_rate units="1/min">0.1</decay_rate>
			</physical_parameter_set>
			<initial_condition units="mmHg">38.0</initial_condition>
			<Dirichlet_boundary_condition units="mmHg" enabled="True">21.0</Dirichlet_boundary_condition>
		</variable>
		<options>
			<calculate_gradients>false</calculate_gradients>
			<track_internalized_substrates_in_each_agent>false</track_internalized_substrates_in_each_agent>
		</options>
	</microenvironment_setup>
</PhysiCell_settings>
)";
	ofs.close();

	physicell_config config = parse_physicell_config(temp_file);

	// Verify substrate parsed correctly
	ASSERT_EQ(config.microenvironment.variables.size(), 1);
	const auto& oxygen = config.microenvironment.variables[0];
	EXPECT_EQ(oxygen.name, "oxygen");
	EXPECT_EQ(oxygen.diffusion_coefficient, 100000.0);

	// Legacy Dirichlet should apply same value to all boundaries with all enabled
	EXPECT_EQ(oxygen.boundary_conditions.mins_values[0], 21.0);
	EXPECT_EQ(oxygen.boundary_conditions.mins_values[1], 21.0);
	EXPECT_EQ(oxygen.boundary_conditions.mins_values[2], 21.0);
	EXPECT_EQ(oxygen.boundary_conditions.maxs_values[0], 21.0);
	EXPECT_EQ(oxygen.boundary_conditions.maxs_values[1], 21.0);
	EXPECT_EQ(oxygen.boundary_conditions.maxs_values[2], 21.0);

	EXPECT_TRUE(oxygen.boundary_conditions.mins_conditions[0]);
	EXPECT_TRUE(oxygen.boundary_conditions.mins_conditions[1]);
	EXPECT_TRUE(oxygen.boundary_conditions.mins_conditions[2]);
	EXPECT_TRUE(oxygen.boundary_conditions.maxs_conditions[0]);
	EXPECT_TRUE(oxygen.boundary_conditions.maxs_conditions[1]);
	EXPECT_TRUE(oxygen.boundary_conditions.maxs_conditions[2]);

	std::filesystem::remove(temp_file);
}

TEST(ConfigReaderErrorTest, NoDirichletOptionsUsesDefaults)
{
	// Create XML with no Dirichlet options at all
	const std::filesystem::path temp_file = "no_dirichlet_test.xml";
	std::ofstream ofs(temp_file);
	ofs << R"(<?xml version="1.0"?>
<PhysiCell_settings>
	<domain>
		<x_min>-100</x_min>
		<x_max>100</x_max>
		<y_min>-100</y_min>
		<y_max>100</y_max>
		<z_min>-100</z_min>
		<z_max>100</z_max>
		<dx>10</dx>
		<dy>10</dy>
		<dz>10</dz>
		<use_2D>false</use_2D>
	</domain>

	<overall>
		<max_time units="min">100</max_time>
		<time_units>min</time_units>
		<space_units>micron</space_units>
		<dt_diffusion units="min">0.01</dt_diffusion>
		<dt_mechanics units="min">0.1</dt_mechanics>
		<dt_phenotype units="min">6</dt_phenotype>
	</overall>

	<microenvironment_setup>
		<variable name="glucose" units="mM" ID="0">
			<physical_parameter_set>
				<diffusion_coefficient units="micron^2/min">600.0</diffusion_coefficient>
				<decay_rate units="1/min">0.05</decay_rate>
			</physical_parameter_set>
			<initial_condition units="mM">5.0</initial_condition>
		</variable>
		<options>
			<calculate_gradients>true</calculate_gradients>
			<track_internalized_substrates_in_each_agent>true</track_internalized_substrates_in_each_agent>
		</options>
	</microenvironment_setup>
</PhysiCell_settings>
)";
	ofs.close();

	physicell_config config = parse_physicell_config(temp_file);

	// Verify substrate parsed correctly
	ASSERT_EQ(config.microenvironment.variables.size(), 1);
	const auto& glucose = config.microenvironment.variables[0];
	EXPECT_EQ(glucose.name, "glucose");

	// Should have default values (all zeros and disabled)
	EXPECT_EQ(glucose.boundary_conditions.mins_values[0], 0.0);
	EXPECT_EQ(glucose.boundary_conditions.mins_values[1], 0.0);
	EXPECT_EQ(glucose.boundary_conditions.mins_values[2], 0.0);
	EXPECT_EQ(glucose.boundary_conditions.maxs_values[0], 0.0);
	EXPECT_EQ(glucose.boundary_conditions.maxs_values[1], 0.0);
	EXPECT_EQ(glucose.boundary_conditions.maxs_values[2], 0.0);

	EXPECT_FALSE(glucose.boundary_conditions.mins_conditions[0]);
	EXPECT_FALSE(glucose.boundary_conditions.mins_conditions[1]);
	EXPECT_FALSE(glucose.boundary_conditions.mins_conditions[2]);
	EXPECT_FALSE(glucose.boundary_conditions.maxs_conditions[0]);
	EXPECT_FALSE(glucose.boundary_conditions.maxs_conditions[1]);
	EXPECT_FALSE(glucose.boundary_conditions.maxs_conditions[2]);

	std::filesystem::remove(temp_file);
}

TEST(ConfigReaderSolverTest, SolverConfigParsed)
{
	// Create XML with solver section
	const std::filesystem::path temp_file = "solver_test.xml";
	std::ofstream ofs(temp_file);
	ofs << R"(<?xml version="1.0"?>
<PhysiCell_settings>
	<domain>
		<x_min>-100</x_min>
		<x_max>100</x_max>
		<y_min>-100</y_min>
		<y_max>100</y_max>
		<z_min>-100</z_min>
		<z_max>100</z_max>
		<dx>10</dx>
		<dy>10</dy>
		<dz>10</dz>
		<use_2D>false</use_2D>
	</domain>

	<overall>
		<max_time units="min">100</max_time>
		<time_units>min</time_units>
		<space_units>micron</space_units>
		<dt_diffusion units="min">0.01</dt_diffusion>
		<dt_mechanics units="min">0.1</dt_mechanics>
		<dt_phenotype units="min">6</dt_phenotype>
	</overall>

	<microenvironment_setup>
		<variable name="oxygen" units="mmHg" ID="0">
			<physical_parameter_set>
				<diffusion_coefficient units="micron^2/min">100000.0</diffusion_coefficient>
				<decay_rate units="1/min">0.1</decay_rate>
			</physical_parameter_set>
			<initial_condition units="mmHg">38.0</initial_condition>
		</variable>
		<options>
			<calculate_gradients>false</calculate_gradients>
			<track_internalized_substrates_in_each_agent>false</track_internalized_substrates_in_each_agent>
		</options>
	</microenvironment_setup>

	<solver>
		<name>openmp</name>
	</solver>
</PhysiCell_settings>
)";
	ofs.close();

	const physicell_config config = parse_physicell_config(temp_file);

	// Verify solver parsed correctly
	EXPECT_EQ(config.solver.name, "openmp");

	std::filesystem::remove(temp_file);
}

TEST(ConfigReaderSolverTest, MissingSolverUsesDefault)
{
	// Create XML without solver section
	const std::filesystem::path temp_file = "no_solver_test.xml";
	std::ofstream ofs(temp_file);
	ofs << R"(<?xml version="1.0"?>
<PhysiCell_settings>
	<domain>
		<x_min>-100</x_min>
		<x_max>100</x_max>
		<y_min>-100</y_min>
		<y_max>100</y_max>
		<z_min>-100</z_min>
		<z_max>100</z_max>
		<dx>10</dx>
		<dy>10</dy>
		<dz>10</dz>
		<use_2D>false</use_2D>
	</domain>

	<overall>
		<max_time units="min">100</max_time>
		<time_units>min</time_units>
		<space_units>micron</space_units>
		<dt_diffusion units="min">0.01</dt_diffusion>
		<dt_mechanics units="min">0.1</dt_mechanics>
		<dt_phenotype units="min">6</dt_phenotype>
	</overall>

	<microenvironment_setup>
		<variable name="oxygen" units="mmHg" ID="0">
			<physical_parameter_set>
				<diffusion_coefficient units="micron^2/min">100000.0</diffusion_coefficient>
				<decay_rate units="1/min">0.1</decay_rate>
			</physical_parameter_set>
			<initial_condition units="mmHg">38.0</initial_condition>
		</variable>
		<options>
			<calculate_gradients>false</calculate_gradients>
			<track_internalized_substrates_in_each_agent>false</track_internalized_substrates_in_each_agent>
		</options>
	</microenvironment_setup>
</PhysiCell_settings>
)";
	ofs.close();

	const physicell_config config = parse_physicell_config(temp_file);

	// Should default to empty string
	EXPECT_EQ(config.solver.name, "");

	std::filesystem::remove(temp_file);
}
