#include <cmath>
#include <memory>

#include <common/base_agent_data.h>
#include <gtest/gtest.h>

#include <micromechanics/agent_container.h>
#include <micromechanics/agent_data.h>
#include <micromechanics/environment.h>
#include <micromechanics/simulation_parameters.h>

#include "potentials/kelvin_voigt_potential.h"
#include "potentials/morse_potential.h"

using namespace physicore;
using namespace physicore::mechanics::micromechanics;
using namespace physicore::mechanics::micromechanics::kernels::openmp_solver;

class PotentialTest : public ::testing::Test
{
protected:
	static std::unique_ptr<environment> create_env_with_two_agents(real_t pos_i_x, real_t pos_j_x)
	{
		auto env = std::make_unique<environment>(0.01);
		auto base_data = std::make_unique<base_agent_data>(3);
		auto mech_data = std::make_unique<agent_data>(*base_data);
		env->agents = std::make_unique<agent_container>(std::move(base_data), std::move(mech_data));

		// Create two agents
		env->agents->create();
		env->agents->create();

		auto& agents = *env->agents;
		auto& base = *std::get<std::unique_ptr<base_agent_data>>(agents.agent_datas);
		auto& mech = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);

		// Agent 0 at (pos_i_x, 0, 0)
		base.positions[0] = pos_i_x;
		base.positions[1] = 0.0;
		base.positions[2] = 0.0;

		// Agent 1 at (pos_j_x, 0, 0)
		base.positions[3] = pos_j_x;
		base.positions[4] = 0.0;
		base.positions[5] = 0.0;

		// Set cell_ids to valid cells
		mech.cell_ids[0] = 0;
		mech.cell_ids[1] = 0;

		// Set up cell data with one cell
		env->cells.resize(1, 3);
		env->cells.radii[0] = 1.0;

		return env;
	}
};

// === Morse Potential Tests ===

TEST_F(PotentialTest, MorseForceZeroAtEquilibrium)
{
	interaction_config config;
	config.morse_scaling_factor = 1.0;
	config.morse_equilibrium_distance = 2.0;
	config.morse_stiffness = 1.0;

	morse_potential pot(config);

	auto env = create_env_with_two_agents(0.0, 2.0);
	// At r = r0, exp_power = a*(1 - r²/r₀²) = 1*(1 - 4/4) = 0
	// exp(0) = 1, so (exp(2*0) - exp(0)) = (1 - 1) = 0 → force = 0
	real_t force = pot.calculate_pairwise_force(*env, 0, 1, 2.0, 2.0, 0.0, 0.0);
	EXPECT_NEAR(force, 0.0, 1e-12);
}

TEST_F(PotentialTest, MorseForceRepulsiveWhenClose)
{
	interaction_config config;
	config.morse_scaling_factor = 1.0;
	config.morse_equilibrium_distance = 2.0;
	config.morse_stiffness = 1.0;

	morse_potential pot(config);

	auto env = create_env_with_two_agents(0.0, 0.5);
	// r=0.5 < r0=2.0, so r²/r₀² < 1, P > 0, exp(2P) > exp(P) → force > 0 (repulsive)
	real_t force = pot.calculate_pairwise_force(*env, 0, 1, 0.5, 0.5, 0.0, 0.0);
	EXPECT_GT(force, 0.0);
}

TEST_F(PotentialTest, MorseForceAttractiveWhenFar)
{
	interaction_config config;
	config.morse_scaling_factor = 1.0;
	config.morse_equilibrium_distance = 2.0;
	config.morse_stiffness = 1.0;

	morse_potential pot(config);

	auto env = create_env_with_two_agents(0.0, 3.0);
	// r=3.0 > r0=2.0, so r²/r₀² > 1, P < 0, exp(2P) < exp(P) → force < 0 (attractive)
	real_t force = pot.calculate_pairwise_force(*env, 0, 1, 3.0, 3.0, 0.0, 0.0);
	EXPECT_LT(force, 0.0);
}

TEST_F(PotentialTest, MorseForceZeroWhenParamsZero)
{
	interaction_config config;
	config.morse_scaling_factor = 0.0;
	config.morse_equilibrium_distance = 2.0;
	config.morse_stiffness = 1.0;

	morse_potential pot(config);

	auto env = create_env_with_two_agents(0.0, 1.0);
	real_t force = pot.calculate_pairwise_force(*env, 0, 1, 1.0, 1.0, 0.0, 0.0);
	EXPECT_DOUBLE_EQ(force, 0.0);
}

TEST_F(PotentialTest, MorseMaxInteractionDistance)
{
	interaction_config config;
	config.morse_equilibrium_distance = 2.0;

	morse_potential pot(config);
	auto env = create_env_with_two_agents(0.0, 1.0);
	EXPECT_DOUBLE_EQ(pot.max_interaction_distance(*env, 0), 5.0); // 2.0 * 2.5
}

// === Kelvin-Voigt Potential Tests ===

TEST_F(PotentialTest, KelvinVoigtSpringForceAtRest)
{
	interaction_config config;
	config.spring_constant = 1.0;
	config.damping_coefficient = 0.0; // no damping

	kelvin_voigt_potential pot(config);

	// radius=1.0, rest_length = 2*radius = 2.0
	auto env = create_env_with_two_agents(0.0, 2.0);
	// At rest length, spring force = k * (distance - rest) = 1.0 * (2.0 - 2.0) = 0
	real_t force = pot.calculate_pairwise_force(*env, 0, 1, 2.0, 2.0, 0.0, 0.0);
	EXPECT_NEAR(force, 0.0, 1e-12);
}

TEST_F(PotentialTest, KelvinVoigtSpringForceCompression)
{
	interaction_config config;
	config.spring_constant = 2.0;
	config.damping_coefficient = 0.0;

	kelvin_voigt_potential pot(config);

	// radius=1.0, rest_length=2.0, distance=1.0 → F = 2.0*(1.0 - 2.0) = -2.0
	auto env = create_env_with_two_agents(0.0, 1.0);
	real_t force = pot.calculate_pairwise_force(*env, 0, 1, 1.0, 1.0, 0.0, 0.0);
	EXPECT_DOUBLE_EQ(force, -2.0);
}

