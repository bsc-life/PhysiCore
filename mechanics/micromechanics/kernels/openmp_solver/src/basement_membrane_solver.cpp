#include "basement_membrane_solver.h"

#include <algorithm>
#include <cmath>
#include <tuple>

#include <micromechanics/agent_data.h>
#include <micromechanics/environment.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void basement_membrane_solver::initialize(environment& /*e*/)
{
	if (initialized_)
		return;

	initialized_ = true;
}

void basement_membrane_solver::update_interactions(environment& e)
{
	if (!e.params.enable_basement_membrane)
		return;

	auto& agents = *e.agents;
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);
	auto& base_data = mech_data.base_data;
	index_t count = agents.size();

	real_t repulsion_strength = e.params.cell_BM_repulsion_strength;

	// Domain boundaries
	real_t x_min = e.domain_min[0];
	real_t x_max = e.domain_max[0];
	real_t y_min = e.domain_min[1];
	real_t y_max = e.domain_max[1];
	real_t z_min = e.domain_min[2];
	real_t z_max = e.domain_max[2];

#pragma omp parallel for
	for (index_t i = 0; i < count; ++i)
	{
		if (!mech_data.is_movable[i])
			continue;

		real_t radius = mech_data.radii[i];
		real_t x = base_data.positions[i * 3];
		real_t y = base_data.positions[i * 3 + 1];
		real_t z = base_data.positions[i * 3 + 2];

		// Check each boundary and apply repulsive force if within radius
		// Force = repulsion_strength * (1 - distance/radius)^2

		// X boundaries
		real_t dist_to_x_min = x - x_min;
		if (dist_to_x_min < radius)
		{
			real_t overlap = 1.0 - dist_to_x_min / radius;
			mech_data.forces[i * 3] += repulsion_strength * overlap * overlap;
		}

		real_t dist_to_x_max = x_max - x;
		if (dist_to_x_max < radius)
		{
			real_t overlap = 1.0 - dist_to_x_max / radius;
			mech_data.forces[i * 3] -= repulsion_strength * overlap * overlap;
		}

		// Y boundaries
		real_t dist_to_y_min = y - y_min;
		if (dist_to_y_min < radius)
		{
			real_t overlap = 1.0 - dist_to_y_min / radius;
			mech_data.forces[i * 3 + 1] += repulsion_strength * overlap * overlap;
		}

		real_t dist_to_y_max = y_max - y;
		if (dist_to_y_max < radius)
		{
			real_t overlap = 1.0 - dist_to_y_max / radius;
			mech_data.forces[i * 3 + 1] -= repulsion_strength * overlap * overlap;
		}

		// Z boundaries (only in 3D)
		if (e.params.dims == 3)
		{
			real_t dist_to_z_min = z - z_min;
			if (dist_to_z_min < radius)
			{
				real_t overlap = 1.0 - dist_to_z_min / radius;
				mech_data.forces[i * 3 + 2] += repulsion_strength * overlap * overlap;
			}

			real_t dist_to_z_max = z_max - z;
			if (dist_to_z_max < radius)
			{
				real_t overlap = 1.0 - dist_to_z_max / radius;
				mech_data.forces[i * 3 + 2] -= repulsion_strength * overlap * overlap;
			}
		}
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver

