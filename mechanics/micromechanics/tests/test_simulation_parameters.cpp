#include <cstdint>
#include <utility>

#include <gtest/gtest.h>

#include "micromechanics/simulation_parameters.h"

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

TEST(SimulationParametersTest, AddInteractionSymmetric)
{
	simulation_parameters params;
	interaction_config config;
	config.adhesion_strength = 2.5;

	params.add_interaction(1, 2, config);

	// Both (1,2) and (2,1) should be stored
	const std::pair<std::uint8_t, std::uint8_t> key_12 { 1, 2 };
	const std::pair<std::uint8_t, std::uint8_t> key_21 { 2, 1 };
	ASSERT_NE(params.interactions.find(key_12), params.interactions.end());
	ASSERT_NE(params.interactions.find(key_21), params.interactions.end());
	EXPECT_DOUBLE_EQ(params.interactions[key_12].adhesion_strength, 2.5);
	EXPECT_DOUBLE_EQ(params.interactions[key_21].adhesion_strength, 2.5);
}

TEST(SimulationParametersTest, AddInteractionSelfPair)
{
	simulation_parameters params;
	interaction_config config;
	config.repulsion_strength = 7.0;

	params.add_interaction(3, 3, config);

	// Only one entry for (3,3), no duplicate
	const std::pair<std::uint8_t, std::uint8_t> key_33 { 3, 3 };
	ASSERT_NE(params.interactions.find(key_33), params.interactions.end());
	EXPECT_DOUBLE_EQ(params.interactions[key_33].repulsion_strength, 7.0);
	EXPECT_EQ(params.interactions.size(), 1U);
}

TEST(SimulationParametersTest, GetInteractionFound)
{
	simulation_parameters params;
	interaction_config config;
	config.spring_constant = 42.0;

	params.add_interaction(0, 1, config);

	const auto& retrieved = params.get_interaction(0, 1);
	EXPECT_DOUBLE_EQ(retrieved.spring_constant, 42.0);

	// Symmetric lookup
	const auto& retrieved_sym = params.get_interaction(1, 0);
	EXPECT_DOUBLE_EQ(retrieved_sym.spring_constant, 42.0);
}

TEST(SimulationParametersTest, GetInteractionFallsBackToDefault)
{
	simulation_parameters params;
	params.default_interaction.morse_scaling_factor = 99.0;

	// No interaction registered for (5,6), should return default
	const auto& retrieved = params.get_interaction(5, 6);
	EXPECT_DOUBLE_EQ(retrieved.morse_scaling_factor, 99.0);
}

TEST(SimulationParametersTest, SetSingleTypeInteraction)
{
	simulation_parameters params;
	interaction_config config;
	config.damping_coefficient = 3.14;

	params.set_single_type_interaction(config);

	const std::pair<std::uint8_t, std::uint8_t> key_00 { 0, 0 };
	ASSERT_NE(params.interactions.find(key_00), params.interactions.end());
	EXPECT_DOUBLE_EQ(params.interactions[key_00].damping_coefficient, 3.14);
}

TEST(SimulationParametersTest, DefaultValues)
{
	const simulation_parameters params;

	EXPECT_EQ(params.solver_name, "openmp_solver");
	EXPECT_DOUBLE_EQ(params.mechanics_timestep, 0.1);
	EXPECT_EQ(params.dims, 3);
	EXPECT_TRUE(params.enable_motility);
	EXPECT_FALSE(params.enable_basement_membrane);
	EXPECT_FALSE(params.enable_spring_attachments);
	EXPECT_DOUBLE_EQ(params.cell_BM_repulsion_strength, 10.0);
	EXPECT_EQ(params.maximum_number_of_attachments, 12);
	EXPECT_DOUBLE_EQ(params.attachment_elastic_constant, 0.01);
	EXPECT_DOUBLE_EQ(params.attachment_rate, 0.0);
	EXPECT_DOUBLE_EQ(params.detachment_rate, 0.0);
}

TEST(SimulationParametersTest, InteractionConfigDefaults)
{
	const interaction_config config;

	EXPECT_EQ(config.potential_name, "morse");
	EXPECT_DOUBLE_EQ(config.adhesion_strength, 0.4);
	EXPECT_DOUBLE_EQ(config.repulsion_strength, 10.0);
	EXPECT_DOUBLE_EQ(config.relative_maximum_adhesion_distance, 1.25);
	EXPECT_DOUBLE_EQ(config.spring_constant, 1.0);
	EXPECT_DOUBLE_EQ(config.damping_coefficient, 0.1);
	EXPECT_DOUBLE_EQ(config.morse_scaling_factor, 1.0);
	EXPECT_DOUBLE_EQ(config.morse_equilibrium_distance, 1.0);
	EXPECT_DOUBLE_EQ(config.morse_stiffness, 1.0);
}
