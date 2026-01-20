#include "position_solver.h"

#include <algorithm>
#include <cmath>
#include <tuple>

#include <micromechanics/agent_container.h>
#include <micromechanics/agent_data.h>
#include <micromechanics/environment.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void position_solver::initialize(environment& /*e*/)
{
	if (initialized_)
		return;

	initialized_ = true;
}

void position_solver::update_positions(environment& e)
{
	auto& agents = *e.agents;
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);
	auto& base_data = mech_data.base_data;
	index_t const count = agents.size();
	real_t const dt = e.timestep;

	// Adams-Bashforth 2nd order integration:
	// x_new = x_old + dt * (1.5 * v_new - 0.5 * v_old)
	//
	// This provides better stability than simple Euler integration
	// by using a weighted average of current and previous velocities.

#pragma omp parallel for
	for (index_t i = 0; i < count; ++i)
	{
		if (!mech_data.is_movable[i])
			continue;

		// Current velocity = force / drag (assuming unit mass, overdamped)
		// In PhysiCell, velocity = force / drag_coefficient
		// For simplicity, we use force directly as velocity (drag = 1)
		real_t const vx = mech_data.forces[i * 3];
		real_t const vy = mech_data.forces[i * 3 + 1];
		real_t const vz = mech_data.forces[i * 3 + 2];

		// Get previous velocity
		real_t const vx_old = mech_data.previous_velocities[i * 3];
		real_t const vy_old = mech_data.previous_velocities[i * 3 + 1];
		real_t const vz_old = mech_data.previous_velocities[i * 3 + 2];

		// Adams-Bashforth coefficients
		constexpr real_t ab_new = 1.5;
		constexpr real_t ab_old = 0.5;

		// Calculate displacement
		real_t const dx = dt * (ab_new * vx - ab_old * vx_old);
		real_t const dy = dt * (ab_new * vy - ab_old * vy_old);
		real_t const dz = dt * (ab_new * vz - ab_old * vz_old);

		// Update position
		base_data.positions[i * 3] += dx;
		base_data.positions[i * 3 + 1] += dy;
		base_data.positions[i * 3 + 2] += dz;

		// Store current velocity as previous for next timestep
		mech_data.previous_velocities[i * 3] = vx;
		mech_data.previous_velocities[i * 3 + 1] = vy;
		mech_data.previous_velocities[i * 3 + 2] = vz;
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
