#include <cmath>
#include <memory>

#include <common/base_agent_data.h>
#include <gtest/gtest.h>

#include "micromechanics/agent_container.h"
#include "micromechanics/agent_data.h"
#include "micromechanics/environment.h"
#include "micromechanics/potential_interface.h"
#include "micromechanics/simulation_parameters.h"

// Include internal potential implementations for testing
#include "potentials/kelvin_voigt_potential.h"
#include "potentials/morse_potential.h"
#include "potentials/standard_potential.h"

using namespace physicore;
using namespace physicore::mechanics::micromechanics;
using namespace physicore::mechanics::micromechanics::kernels::openmp_solver;

class PotentialTest : public ::testing::Test
{
protected:
	void SetUp() override
	{
		env = std::make_unique<environment>(0.01);
		auto base_data = std::make_unique<base_agent_data>(3);
		auto mech_data = std::make_unique<agent_data>(*base_data);
		env->agents = std::make_unique<agent_container>(std::move(base_data), std::move(mech_data));
	}

	void AddAgent(real_t x, real_t y, real_t z, real_t radius)
	{
		auto* agent = env->agents->create();
		agent->position()[0] = x;
		agent->position()[1] = y;
		agent->position()[2] = z;
		agent->radius() = radius;
	}

	std::unique_ptr<environment> env;
};

// ============== Standard Potential Tests ==============

TEST_F(PotentialTest, StandardPotentialRepulsion)
{
	// Setup: two overlapping agents
	AddAgent(0, 0, 0, 10);
	AddAgent(15, 0, 0, 10); // Distance 15, sum of radii = 20 -> overlap

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
	mech_data.cell_cell_repulsion_strength[0] = 10.0;
	mech_data.cell_cell_repulsion_strength[1] = 10.0;
	mech_data.cell_cell_adhesion_strength[0] = 0.0;
	mech_data.cell_cell_adhesion_strength[1] = 0.0;

	interaction_config const config;
	standard_potential const potential(config);

	real_t force_out = 0.0;
	real_t const distance = 15.0;
	potential.calculate_pairwise_force(*env, 0, 1, distance, 15.0, 0.0, 0.0, force_out);

	// Repulsion should be positive (pushing apart)
	EXPECT_GT(force_out, 0.0);
}

TEST_F(PotentialTest, StandardPotentialAdhesion)
{
	// Setup: two agents just touching (no overlap, but within adhesion range)
	AddAgent(0, 0, 0, 10);
	AddAgent(22, 0, 0, 10); // Distance 22, sum of radii = 20 -> no overlap, slight separation

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
	mech_data.cell_cell_repulsion_strength[0] = 10.0;
	mech_data.cell_cell_repulsion_strength[1] = 10.0;
	mech_data.cell_cell_adhesion_strength[0] = 10.0;
	mech_data.cell_cell_adhesion_strength[1] = 10.0;
	mech_data.relative_maximum_adhesion_distance[0] = 1.5;
	mech_data.relative_maximum_adhesion_distance[1] = 1.5;

	interaction_config const config;
	standard_potential const potential(config);

	real_t force_out = 0.0;
	real_t const distance = 22.0;
	potential.calculate_pairwise_force(*env, 0, 1, distance, 22.0, 0.0, 0.0, force_out);

	// Adhesion should dominate (negative force = pulling together)
	EXPECT_LT(force_out, 0.0);
}

TEST_F(PotentialTest, StandardPotentialNoForceOutOfRange)
{
	// Setup: two agents far apart
	AddAgent(0, 0, 0, 10);
	AddAgent(100, 0, 0, 10); // Distance 100, far out of range

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
	mech_data.cell_cell_repulsion_strength[0] = 10.0;
	mech_data.cell_cell_repulsion_strength[1] = 10.0;
	mech_data.cell_cell_adhesion_strength[0] = 10.0;
	mech_data.cell_cell_adhesion_strength[1] = 10.0;
	mech_data.relative_maximum_adhesion_distance[0] = 1.5; // Max dist = 1.5*10 + 1.5*10 = 30
	mech_data.relative_maximum_adhesion_distance[1] = 1.5;

	interaction_config const config;
	standard_potential const potential(config);

	real_t force_out = 0.0;
	real_t const distance = 100.0;
	potential.calculate_pairwise_force(*env, 0, 1, distance, 100.0, 0.0, 0.0, force_out);

	// No force at this distance
	EXPECT_DOUBLE_EQ(force_out, 0.0);
}

// ============== Morse Potential Tests ==============

TEST_F(PotentialTest, MorsePotentialAtEquilibrium)
{
	AddAgent(0, 0, 0, 10);
	AddAgent(20, 0, 0, 10);

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
	mech_data.intra_scaling_factors[0] = 1.0;
	mech_data.intra_equilibrium_distances[0] = 20.0; // Equilibrium at distance 20
	mech_data.intra_stiffnesses[0] = 1.0;

	interaction_config const config;
	morse_potential const potential(config);

	real_t force_out = 0.0;
	real_t const distance = 20.0;
	potential.calculate_pairwise_force(*env, 0, 1, distance, 20.0, 0.0, 0.0, force_out);

	// At equilibrium, force should be zero
	EXPECT_NEAR(force_out, 0.0, 1e-10);
}

TEST_F(PotentialTest, MorsePotentialRepulsionWhenCompressed)
{
	AddAgent(0, 0, 0, 10);
	AddAgent(15, 0, 0, 10);

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
	mech_data.intra_scaling_factors[0] = 1.0;
	mech_data.intra_equilibrium_distances[0] = 20.0;
	mech_data.intra_stiffnesses[0] = 1.0;

	interaction_config const config;
	morse_potential const potential(config);

	real_t force_out = 0.0;
	real_t const distance = 15.0; // Less than equilibrium
	potential.calculate_pairwise_force(*env, 0, 1, distance, 15.0, 0.0, 0.0, force_out);

	// When compressed (r < r0), should repel (positive force)
	EXPECT_GT(force_out, 0.0);
}

TEST_F(PotentialTest, MorsePotentialAttractionWhenStretched)
{
	AddAgent(0, 0, 0, 10);
	AddAgent(25, 0, 0, 10);

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
	mech_data.intra_scaling_factors[0] = 1.0;
	mech_data.intra_equilibrium_distances[0] = 20.0;
	mech_data.intra_stiffnesses[0] = 1.0;

	interaction_config const config;
	morse_potential const potential(config);

	real_t force_out = 0.0;
	real_t const distance = 25.0; // Greater than equilibrium
	potential.calculate_pairwise_force(*env, 0, 1, distance, 25.0, 0.0, 0.0, force_out);

	// When stretched (r > r0), should attract (negative force)
	EXPECT_LT(force_out, 0.0);
}

TEST_F(PotentialTest, StandardPotentialBalancePoint)
{
	// Test where repulsion and adhesion balance
	AddAgent(0, 0, 0, 10);
	AddAgent(20, 0, 0, 10); // Exactly at sum of radii (touching)

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
	// Equal repulsion and adhesion strengths
	mech_data.cell_cell_repulsion_strength[0] = 10.0;
	mech_data.cell_cell_repulsion_strength[1] = 10.0;
	mech_data.cell_cell_adhesion_strength[0] = 10.0;
	mech_data.cell_cell_adhesion_strength[1] = 10.0;
	mech_data.relative_maximum_adhesion_distance[0] = 1.5;
	mech_data.relative_maximum_adhesion_distance[1] = 1.5;

	interaction_config const config;
	standard_potential const potential(config);

	real_t force_out = 0.0;
	real_t const distance = 20.0; // Exactly at repulsive boundary
	potential.calculate_pairwise_force(*env, 0, 1, distance, 20.0, 0.0, 0.0, force_out);

	// At touching point: repulsion=0 (1-d/R = 0), but adhesion > 0
	// So net force should be negative (attraction)
	EXPECT_LT(force_out, 0.0);
}

// ============== Kelvin-Voigt Potential Tests ==============

TEST_F(PotentialTest, KelvinVoigtSpringForceAtRest)
{
	AddAgent(0, 0, 0, 10);
	AddAgent(20, 0, 0, 10); // Rest length = 2*radius = 20

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
	mech_data.spring_constants[0] = 5.0;
	mech_data.dissipation_rates[0] = 0.0; // No damping for this test
	mech_data.radii[0] = 10.0;

	interaction_config const config;
	kelvin_voigt_potential const potential(config);

	real_t force_out = 0.0;
	real_t const distance = 20.0; // At rest length
	potential.calculate_pairwise_force(*env, 0, 1, distance, 20.0, 0.0, 0.0, force_out);

	// At rest length, spring force should be zero
	EXPECT_NEAR(force_out, 0.0, 1e-10);
}

TEST_F(PotentialTest, KelvinVoigtSpringForceCompressed)
{
	AddAgent(0, 0, 0, 10);
	AddAgent(15, 0, 0, 10);

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
	mech_data.spring_constants[0] = 5.0;
	mech_data.dissipation_rates[0] = 0.0; // No damping
	mech_data.radii[0] = 10.0;			  // Rest length = 20

	interaction_config const config;
	kelvin_voigt_potential const potential(config);

	real_t force_out = 0.0;
	real_t const distance = 15.0; // Compressed by 5 units
	potential.calculate_pairwise_force(*env, 0, 1, distance, 15.0, 0.0, 0.0, force_out);

	// F = k * (d - rest) = 5 * (15 - 20) = -25 (negative = push apart)
	EXPECT_DOUBLE_EQ(force_out, -25.0);
}

TEST_F(PotentialTest, KelvinVoigtSpringForceStretched)
{
	AddAgent(0, 0, 0, 10);
	AddAgent(25, 0, 0, 10);

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
	mech_data.spring_constants[0] = 5.0;
	mech_data.dissipation_rates[0] = 0.0; // No damping
	mech_data.radii[0] = 10.0;			  // Rest length = 20

	interaction_config const config;
	kelvin_voigt_potential const potential(config);

	real_t force_out = 0.0;
	real_t const distance = 25.0; // Stretched by 5 units
	potential.calculate_pairwise_force(*env, 0, 1, distance, 25.0, 0.0, 0.0, force_out);

	// F = k * (d - rest) = 5 * (25 - 20) = 25 (positive = pull together)
	EXPECT_DOUBLE_EQ(force_out, 25.0);
}

TEST_F(PotentialTest, KelvinVoigtDampingForce)
{
	AddAgent(0, 0, 0, 10);
	AddAgent(20, 0, 0, 10);

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
	mech_data.spring_constants[0] = 0.0; // No spring for this test
	mech_data.dissipation_rates[0] = 2.0;
	mech_data.radii[0] = 10.0;

	// Set velocities: agent 1 moving toward agent 0
	mech_data.previous_velocities[0] = 0.0;	 // agent 0 stationary
	mech_data.previous_velocities[3] = -5.0; // agent 1 moving left (toward agent 0)

	interaction_config const config;
	kelvin_voigt_potential const potential(config);

	real_t force_out = 0.0;
	real_t const distance = 20.0;
	real_t const dx = 20.0; // Direction from agent 0 to agent 1
	potential.calculate_pairwise_force(*env, 0, 1, distance, dx, 0.0, 0.0, force_out);

	// dv = v1 - v0 = -5 - 0 = -5 (relative velocity toward agent 0)
	// v_rel_dot_n = dv * dx = -5 * 20 = -100
	// F_damp = gamma * dt * v_rel_dot_n = 2 * 0.01 * (-100) = -2
	EXPECT_DOUBLE_EQ(force_out, -2.0);
}

TEST_F(PotentialTest, KelvinVoigtCombinedForce)
{
	AddAgent(0, 0, 0, 10);
	AddAgent(25, 0, 0, 10);

	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(env->agents->agent_datas);
	mech_data.spring_constants[0] = 5.0;
	mech_data.dissipation_rates[0] = 2.0;
	mech_data.radii[0] = 10.0;

	// Set velocities: agents moving apart
	mech_data.previous_velocities[0] = -1.0; // agent 0 moving left
	mech_data.previous_velocities[3] = 1.0;	 // agent 1 moving right

	interaction_config const config;
	kelvin_voigt_potential const potential(config);

	real_t force_out = 0.0;
	real_t const distance = 25.0;
	real_t const dx = 25.0;
	potential.calculate_pairwise_force(*env, 0, 1, distance, dx, 0.0, 0.0, force_out);

	// Spring: F = k * (d - rest) = 5 * (25 - 20) = 25
	// Damping: dv = 1 - (-1) = 2, v_dot_n = 2 * 25 = 50
	//          F_damp = 2 * 0.01 * 50 = 1
	// Total = 25 + 1 = 26
	EXPECT_DOUBLE_EQ(force_out, 26.0);
}
