#include "kelvin_voigt_potential.h"

#include <utility>

#include <micromechanics/agent_container.h>
#include <micromechanics/agent_data.h>
#include <micromechanics/environment.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

kelvin_voigt_potential::kelvin_voigt_potential(interaction_config config) : config_(std::move(config)) {}

real_t kelvin_voigt_potential::calculate_pairwise_force(const environment& env, index_t agent_i, index_t agent_j,
														real_t distance, real_t dx, real_t dy, real_t dz) const
{
	auto& agents = *env.agents;
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);

	// All spring parameters from config (type-pair level)
	real_t const spring_constant = config_.spring_constant;
	real_t const damping = config_.damping_coefficient;

	// Rest length: look up radii from cell_data via cell_id
	index_t const cell_id_i = mech_data.cell_ids[static_cast<std::size_t>(agent_i)];
	real_t const r_i =
		(cell_id_i != cell_data::invalid_cell_id) ? env.cells.radii[static_cast<std::size_t>(cell_id_i)] : 1.0;
	real_t const rest_length = r_i * 2.0;

	// Spring force: F_spring = k * (distance - rest_length)
	real_t const force_spring = spring_constant * (distance - rest_length);

	// Damping force: F_damp = gamma * (v_rel . n) * dt
	index_t const dims = 3;
	real_t const dvx =
		mech_data.previous_velocities[agent_j * dims + 0] - mech_data.previous_velocities[agent_i * dims + 0];
	real_t const dvy =
		mech_data.previous_velocities[agent_j * dims + 1] - mech_data.previous_velocities[agent_i * dims + 1];
	real_t const dvz =
		mech_data.previous_velocities[agent_j * dims + 2] - mech_data.previous_velocities[agent_i * dims + 2];

	// Project velocity difference onto normal direction
	real_t const v_rel_dot_n = (dvx * dx + dvy * dy + dvz * dz);

	// Include timestep factor as in legacy code
	real_t const dt = env.timestep;
	real_t const force_damp = damping * dt * v_rel_dot_n;

	// Total force (positive = stretching/repulsion, negative = compression/attraction)
	return force_spring + force_damp;
}

std::string kelvin_voigt_potential::name() const { return "kelvin_voigt"; }

real_t kelvin_voigt_potential::max_interaction_distance(const environment& env, index_t agent_i) const
{
	auto& agents = *env.agents;
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);

	// Look up radius from cell_data
	index_t const cell_id = mech_data.cell_ids[static_cast<std::size_t>(agent_i)];
	real_t const radius =
		(cell_id != cell_data::invalid_cell_id) ? env.cells.radii[static_cast<std::size_t>(cell_id)] : 1.0;

	return radius * 2.5;
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
