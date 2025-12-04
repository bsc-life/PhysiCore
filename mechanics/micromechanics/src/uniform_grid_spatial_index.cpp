#include "micromechanics/uniform_grid_spatial_index.h"

#include <cmath>
#include <tuple>

#include "micromechanics/agent_data.h"
#include "micromechanics/environment.h"

namespace physicore::mechanics::micromechanics {

uniform_grid_spatial_index::uniform_grid_spatial_index(real_t cell_size) : cell_size(cell_size) {}

void uniform_grid_spatial_index::build(const environment& env)
{
	grid.clear();
	auto& agents = *env.agents;
	auto& mech_data_ptr = std::get<std::unique_ptr<agent_data>>(agents.agent_datas);
	auto& base_data = mech_data_ptr->base_data;
	index_t count = agents.size();

	for (index_t i = 0; i < count; ++i)
	{
		real_t x = base_data.positions[i * 3];
		real_t y = base_data.positions[i * 3 + 1];
		real_t z = base_data.positions[i * 3 + 2];

		grid_key key { static_cast<int>(std::floor(x / cell_size)), static_cast<int>(std::floor(y / cell_size)),
					   static_cast<int>(std::floor(z / cell_size)) };
		grid[key].push_back(i);
	}
}

std::vector<index_t> uniform_grid_spatial_index::query_neighbors(const environment& env, index_t agent_index,
																 real_t radius) const
{
	std::vector<index_t> neighbors;
	auto& agents = *env.agents;
	auto& mech_data_ptr = std::get<std::unique_ptr<agent_data>>(agents.agent_datas);
	auto& base_data = mech_data_ptr->base_data;

	real_t x = base_data.positions[agent_index * 3];
	real_t y = base_data.positions[agent_index * 3 + 1];
	real_t z = base_data.positions[agent_index * 3 + 2];

	int cx = static_cast<int>(std::floor(x / cell_size));
	int cy = static_cast<int>(std::floor(y / cell_size));
	int cz = static_cast<int>(std::floor(z / cell_size));

	int search_radius = static_cast<int>(std::ceil(radius / cell_size));

	for (int dx = -search_radius; dx <= search_radius; ++dx)
	{
		for (int dy = -search_radius; dy <= search_radius; ++dy)
		{
			for (int dz = -search_radius; dz <= search_radius; ++dz)
			{
				grid_key key { cx + dx, cy + dy, cz + dz };
				auto it = grid.find(key);
				if (it != grid.end())
				{
					for (index_t other_index : it->second)
					{
						if (agent_index == other_index)
							continue;

						real_t ox = base_data.positions[other_index * 3];
						real_t oy = base_data.positions[other_index * 3 + 1];
						real_t oz = base_data.positions[other_index * 3 + 2];

						real_t dist_sq = (x - ox) * (x - ox) + (y - oy) * (y - oy) + (z - oz) * (z - oz);
						if (dist_sq <= radius * radius)
						{
							neighbors.push_back(other_index);
						}
					}
				}
			}
		}
	}
	return neighbors;
}

} // namespace physicore::mechanics::micromechanics
