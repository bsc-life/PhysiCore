#include "physicell/openmp_solver/cell_neighbors.h"

#include <cmath>

namespace physicore::mechanics::physicell::kernels::openmp_solver {

void update_cell_neighbors(environment& e)
{
	auto& data = e.get_agent_data();
	const index_t dims = data.base_data.dims;
	const index_t agents_count = data.agents_count;

	const real_t* positions = data.base_data.positions.data();

	// Clear all neighbor lists first, to match the behavior of the OpenMP solver.
	for (index_t i = 0; i < agents_count; ++i)
		data.state_data.neighbors[i].clear();

	for (index_t i = 0; i < agents_count; ++i)
	{
		if (data.state_data.is_movable[i] == 0)
			continue;

		for (index_t j = 0; j < agents_count; ++j)
		{
			if (j == i)
				continue;

			const real_t adhesion_distance =
				data.mechanics_data.relative_maximum_adhesion_distance[i] * data.radius[i]
				+ data.mechanics_data.relative_maximum_adhesion_distance[j] * data.radius[j];

			if (adhesion_distance <= 0)
				continue;

			real_t distance_sq = 0;
			for (index_t d = 0; d < dims; ++d)
			{
				const real_t diff = positions[i * dims + d] - positions[j * dims + d];
				distance_sq += diff * diff;
			}

			const real_t distance = std::sqrt(distance_sq);
			if (distance <= adhesion_distance)
			{
				data.state_data.neighbors[i].push_back(j);
			}
		}
	}
}

} // namespace physicore::mechanics::physicell::kernels::openmp_solver

