#include "standard_potential.h"

#include <algorithm>
#include <cmath>
#include <tuple>
#include <utility>

#include <micromechanics/agent_data.h>
#include <micromechanics/environment.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

standard_potential::standard_potential(interaction_config config) : config_(std::move(config)) {}

void standard_potential::calculate_pairwise_force(const environment& env, index_t agent_i, index_t agent_j,
												  real_t distance, real_t /*dx*/, real_t /*dy*/, real_t /*dz*/,
												  real_t& force_out) const
{
	auto& agents = *env.agents;
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);

	real_t const radius_i = mech_data.radii[agent_i];
	real_t const radius_j = mech_data.radii[agent_j];

	force_out = 0.0;

	// ========== REPULSION (matches legacy solve_pair exactly) ==========
	// Repulsive distance is sum of radii
	real_t const repulsive_distance = radius_i + radius_j;

	real_t repulsion = 1.0 - distance / repulsive_distance;
	repulsion = std::max(repulsion, 0.0);

	repulsion *= repulsion; // Quadratic falloff

	// Geometric mean of per-agent repulsion strengths (legacy formula)
	real_t const c_rep =
		std::sqrt(mech_data.cell_cell_repulsion_strength[agent_i] * mech_data.cell_cell_repulsion_strength[agent_j]);

	force_out += c_rep * repulsion;

	// ========== ADHESION (matches legacy solve_pair exactly) ==========
	// Adhesion distance uses per-agent relative_maximum_adhesion_distance
	// Legacy: adhesion_dist = rel_max_dist[i] * radius[i] + rel_max_dist[j] * radius[j]
	real_t const adhesion_distance = mech_data.relative_maximum_adhesion_distance[agent_i] * radius_i
									 + mech_data.relative_maximum_adhesion_distance[agent_j] * radius_j;

	real_t adhesion = 1.0 - distance / adhesion_distance;
	adhesion = std::max(adhesion, 0.0);

	adhesion *= adhesion; // Quadratic falloff

	// Geometric mean of per-agent adhesion strengths (legacy formula)
	// Note: Legacy also includes cell_adhesion_affinity matrix - we use config scaling for now
	real_t const c_adh =
		std::sqrt(mech_data.cell_cell_adhesion_strength[agent_i] * mech_data.cell_cell_adhesion_strength[agent_j]);

	force_out -= c_adh * adhesion;

	// Scale by distance to convert to force per unit displacement
	// Legacy: force = (repulsion - adhesion) / distance, then multiplied by position_difference
	// We return the signed magnitude; caller applies direction
	force_out /= distance;
}

std::string standard_potential::name() const { return "standard"; }

real_t standard_potential::max_interaction_distance(const environment& env, index_t agent_i) const
{
	auto& agents = *env.agents;
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);

	// Max distance is when adhesion can still occur
	// Use the agent's own relative_maximum_adhesion_distance * 2 * radius as conservative estimate
	return mech_data.relative_maximum_adhesion_distance[agent_i] * mech_data.radii[agent_i] * 2.0;
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
