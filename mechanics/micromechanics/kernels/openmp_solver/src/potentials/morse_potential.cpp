#include "morse_potential.h"

#include <cmath>
#include <utility>

#include <micromechanics/environment.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

morse_potential::morse_potential(interaction_config config) : config_(std::move(config)) {}

real_t morse_potential::calculate_pairwise_force(const environment& /*env*/, index_t /*agent_i*/, index_t /*agent_j*/,
												 real_t distance, real_t /*dx*/, real_t /*dy*/, real_t /*dz*/) const
{
	// All Morse parameters come from config (type-pair level)
	real_t const scaling_factor = config_.morse_scaling_factor;
	real_t const equilibrium_distance = config_.morse_equilibrium_distance;
	real_t const stiffness = config_.morse_stiffness;

	// Avoid division by zero
	if (scaling_factor == 0.0 || equilibrium_distance == 0.0)
		return 0.0;

	// Potential well depth:  D = (k · r₀²) / (8 · a²)
	real_t const r0_sq = equilibrium_distance * equilibrium_distance;
	real_t const potential_well_depth = (stiffness * r0_sq) / (8.0 * scaling_factor * scaling_factor);

	// Exponent:  P = a · (1 − r²/r₀²)
	real_t const exp_power = scaling_factor * (1.0 - (distance * distance) / r0_sq);

	// Force (−dV/dr):  F = (4·a·r·D / r₀²) · [ exp(2P) − exp(P) ]
	real_t const exp_p = std::exp(exp_power);
	return (4.0 * scaling_factor * distance * potential_well_depth) * (exp_p * exp_p - exp_p) / r0_sq;
}

std::string morse_potential::name() const { return "morse"; }

real_t morse_potential::max_interaction_distance(const environment& /*env*/, index_t /*agent_i*/) const
{
	// Morse potential has longer range — use 2.5× equilibrium distance
	return config_.morse_equilibrium_distance * 2.5;
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
