#include "motility_solver.h"

#include <cmath>
#include <random>
#include <tuple>

#include <micromechanics/agent_container.h>
#include <micromechanics/agent_data.h>
#include <micromechanics/environment.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void motility_solver::initialize(environment& /*e*/)
{
	if (initialized_)
		return;

	initialized_ = true;
}

void motility_solver::update_motility(environment& e)
{
	if (!e.params.enable_motility)
		return;

	auto& agents = *e.agents;
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);
	index_t const count = agents.size();
	real_t const dt = e.timestep;
	auto& cells = e.cells;

	// --- Phase 1: update motility directions per cell ---
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<real_t> uniform(0.0, 1.0);

		for (index_t c = 0; c < cells.cells_count; ++c)
		{
			if (!cells.is_motile[static_cast<std::size_t>(c)])
				continue;

			real_t const persistence_time = cells.persistence_times[static_cast<std::size_t>(c)];
			real_t const migration_bias = cells.migration_biases[static_cast<std::size_t>(c)];

			if (persistence_time > 0.0 && uniform(gen) < dt / persistence_time)
			{
				// Generate random direction on unit sphere
				real_t const theta = 2.0 * M_PI * uniform(gen);
				real_t const phi = std::acos(2.0 * uniform(gen) - 1.0);

				real_t const rand_x = std::sin(phi) * std::cos(theta);
				real_t const rand_y = std::sin(phi) * std::sin(theta);
				real_t const rand_z = std::cos(phi);

				// Get bias direction from cell data
				real_t const bias_x = cells.migration_bias_directions[cells.cell_offset(c, 0)];
				real_t const bias_y = cells.migration_bias_directions[cells.cell_offset(c, 1)];
				real_t const bias_z = cells.migration_bias_directions[cells.cell_offset(c, 2)];

				// Combine random and bias
				real_t dir_x = (1.0 - migration_bias) * rand_x + migration_bias * bias_x;
				real_t dir_y = (1.0 - migration_bias) * rand_y + migration_bias * bias_y;
				real_t dir_z = (1.0 - migration_bias) * rand_z + migration_bias * bias_z;

				// Normalize
				real_t const mag = std::sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
				if (mag > 1e-16)
				{
					dir_x /= mag;
					dir_y /= mag;
					dir_z /= mag;
				}

				cells.motility_directions[cells.cell_offset(c, 0)] = dir_x;
				cells.motility_directions[cells.cell_offset(c, 1)] = dir_y;
				cells.motility_directions[cells.cell_offset(c, 2)] = dir_z;
			}
		}
	}

	// --- Phase 2: apply motility force to each agent ---
#pragma omp parallel for
	for (index_t i = 0; i < count; ++i)
	{
		index_t const cell_id = mech_data.cell_ids[static_cast<std::size_t>(i)];
		if (cell_id == cell_data::invalid_cell_id || cell_id >= cells.cells_count)
			continue;
		if (!cells.is_motile[static_cast<std::size_t>(cell_id)])
			continue;

		real_t const speed = cells.migration_speeds[static_cast<std::size_t>(cell_id)];

		mech_data.forces[i * 3] += speed * cells.motility_directions[cells.cell_offset(cell_id, 0)];
		mech_data.forces[i * 3 + 1] += speed * cells.motility_directions[cells.cell_offset(cell_id, 1)];
		mech_data.forces[i * 3 + 2] += speed * cells.motility_directions[cells.cell_offset(cell_id, 2)];
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
