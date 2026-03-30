#include "micromechanics/cell_aggregation.h"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace physicore::mechanics::micromechanics {

void reset_cell_aggregates(cell_data& cells)
{
	std::fill(cells.positions.begin(), cells.positions.end(), static_cast<real_t>(0.0));
	std::fill(cells.velocities.begin(), cells.velocities.end(), static_cast<real_t>(0.0));
	std::fill(cells.agent_counts.begin(), cells.agent_counts.end(), static_cast<index_t>(0));
	std::fill(cells.pressures.begin(), cells.pressures.end(), static_cast<real_t>(0.0));
}

void aggregate_cell_data_from_agents(const physicore::base_agent_data& base, const agent_data& agents, cell_data& cells)
{
	assert(base.dims == cells.dims);
	assert(agents.base_data.dims == cells.dims);
	const index_t dims = cells.dims;

	reset_cell_aggregates(cells);

	for (index_t a = 0; a < agents.agents_count; ++a)
	{
		const index_t cell_id = agents.cell_ids[static_cast<std::size_t>(a)];
		if (cell_id == cell_data::invalid_cell_id || cell_id >= cells.cells_count)
			continue;

		cells.agent_counts[static_cast<std::size_t>(cell_id)] += 1;

		// COM sums
		for (index_t d = 0; d < dims; ++d)
		{
			cells.positions[cells.cell_offset(cell_id, d)] +=
				base.positions[static_cast<std::size_t>(a) * static_cast<std::size_t>(dims)
							   + static_cast<std::size_t>(d)];
			cells.velocities[cells.cell_offset(cell_id, d)] +=
				agents.velocities[static_cast<std::size_t>(a) * static_cast<std::size_t>(dims)
								  + static_cast<std::size_t>(d)];
		}

		// Pressure proxy: ||F|| (L2 norm of net force on the agent)
		real_t f2 = 0.0;
		for (index_t d = 0; d < dims; ++d)
		{
			const real_t f =
				agents
					.forces[static_cast<std::size_t>(a) * static_cast<std::size_t>(dims) + static_cast<std::size_t>(d)];
			f2 += f * f;
		}
		cells.pressures[static_cast<std::size_t>(cell_id)] += std::sqrt(f2);
	}

	// Turn sums into averages
	for (index_t c = 0; c < cells.cells_count; ++c)
	{
		const index_t n = cells.agent_counts[static_cast<std::size_t>(c)];
		if (n <= 0)
			continue;
		const real_t inv = static_cast<real_t>(1.0) / static_cast<real_t>(n);
		for (index_t d = 0; d < dims; ++d)
		{
			cells.positions[cells.cell_offset(c, d)] *= inv;
			cells.velocities[cells.cell_offset(c, d)] *= inv;
		}
	}
}

} // namespace physicore::mechanics::micromechanics
