#include "morse_potential.h"

#include <cmath>
#include <tuple>
#include <utility>

#include <micromechanics/agent_data.h>
#include <micromechanics/environment.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

morse_potential::morse_potential(interaction_config config) : config_(std::move(config)) {}

void morse_potential::calculate_pairwise_force(const environment& env, index_t agent_i, index_t /*agent_j*/,
											   real_t distance, real_t /*dx*/, real_t /*dy*/, real_t /*dz*/,
											   real_t& force_out) const
{
	auto& agents = *env.agents;
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);

	// Get Morse parameters - prefer per-agent values if set, otherwise use config
	real_t scaling_factor = mech_data.intra_scaling_factors[agent_i];
	real_t equilibrium_distance = mech_data.intra_equilibrium_distances[agent_i];
	real_t stiffness = mech_data.intra_stiffnesses[agent_i];

	// Fall back to config values if agent values are zero
	if (scaling_factor == 0.0)
		scaling_factor = config_.morse_scaling_factor;
	if (equilibrium_distance == 0.0)
		equilibrium_distance = config_.morse_equilibrium_distance;
	if (stiffness == 0.0)
		stiffness = config_.morse_stiffness;

	// Avoid division by zero
	if (scaling_factor == 0.0 || equilibrium_distance == 0.0)
	{
		force_out = 0.0;
		return;
	}

	// Calculate potential well depth
	// D = (k * r0^2) / (8 * a^2)
	real_t const potential_well_depth =
		(stiffness * equilibrium_distance * equilibrium_distance) / (8.0 * scaling_factor * scaling_factor);

	// Calculate exp power: a * (1 - r^2/r0^2)
	real_t const exp_power =
		scaling_factor * (1.0 - (distance * distance) / (equilibrium_distance * equilibrium_distance));

	// Force magnitude from Morse potential derivative
	// F = (4 * a * r * D) * (exp(2P) - exp(P)) / r0^2
	real_t const exp_p = std::exp(exp_power);
	force_out = (4.0 * scaling_factor * distance * potential_well_depth) * (exp_p * exp_p - exp_p)
				/ (equilibrium_distance * equilibrium_distance);

	// Note: positive force_out means repulsion (pushing apart)
	// When r < r0, exp_power > 0, exp(2P) > exp(P), so F > 0 (repulsion)
	// When r > r0, exp_power < 0, exp(2P) < exp(P), so F < 0 (attraction)
}

std::string morse_potential::name() const { return "morse"; }

real_t morse_potential::max_interaction_distance(const environment& env, index_t agent_i) const
{
	auto& agents = *env.agents;
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);

	// Morse potential has longer range - use 2.5x equilibrium distance
	real_t eq_dist = mech_data.intra_equilibrium_distances[agent_i];
	if (eq_dist == 0.0)
		eq_dist = config_.morse_equilibrium_distance;

	return eq_dist * 2.5;
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
