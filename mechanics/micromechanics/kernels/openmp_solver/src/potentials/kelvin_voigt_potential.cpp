#include "kelvin_voigt_potential.h"

#include <cmath>
#include <tuple>
#include <utility>

#include <micromechanics/agent_data.h>
#include <micromechanics/environment.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

kelvin_voigt_potential::kelvin_voigt_potential(interaction_config  config) : config_(std::move(config)) {}

void kelvin_voigt_potential::calculate_pairwise_force(const environment& env, index_t agent_i, index_t agent_j,
													  real_t distance, real_t dx, real_t dy, real_t dz,
													  real_t& force_out) const
{
	auto& agents = *env.agents;
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);

	// Get spring parameters
	real_t spring_constant = mech_data.spring_constants[agent_i];
	if (spring_constant == 0.0)
		spring_constant = config_.spring_constant;

	real_t damping = mech_data.dissipation_rates[agent_i];
	if (damping == 0.0)
		damping = config_.damping_coefficient;

	// Rest length is typically 2 * radius for adjacent cells
	real_t const rest_length = mech_data.radii[agent_i] * 2.0;

	// Spring force: F_spring = k * (distance - rest_length)
	real_t const force_spring = spring_constant * (distance - rest_length);

	// Damping force: F_damp = gamma * (v_rel . n) * dt
	// Get previous velocities
	index_t const dims = 3;
	real_t const dvx = mech_data.previous_velocities[agent_j * dims + 0] - mech_data.previous_velocities[agent_i * dims + 0];
	real_t const dvy = mech_data.previous_velocities[agent_j * dims + 1] - mech_data.previous_velocities[agent_i * dims + 1];
	real_t const dvz = mech_data.previous_velocities[agent_j * dims + 2] - mech_data.previous_velocities[agent_i * dims + 2];

	// Project velocity difference onto normal direction
	real_t const v_rel_dot_n = (dvx * dx + dvy * dy + dvz * dz);

	// Include timestep factor as in legacy code
	real_t const dt = env.timestep;
	real_t const force_damp = damping * dt * v_rel_dot_n;

	// Total force (positive = stretching/repulsion, negative = compression/attraction)
	force_out = force_spring + force_damp;
}

std::string kelvin_voigt_potential::name() const { return "kelvin_voigt"; }

real_t kelvin_voigt_potential::max_interaction_distance(const environment& env, index_t agent_i) const
{
	auto& agents = *env.agents;
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);

	// Kelvin-Voigt typically only applies to connected cells
	// Use 2.5x radius to find potential springs
	return mech_data.radii[agent_i] * 2.5;
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
