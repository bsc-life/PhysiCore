#include "micromechanics/uniform_grid_spatial_index.h"

#include <cmath>
#include <tuple>

#include "micromechanics/agent_container.h"
#include "micromechanics/agent_data.h"
#include "micromechanics/environment.h"

namespace physicore::mechanics::micromechanics {

namespace {

void collect_neighbors_in_cell(const std::vector<index_t>& cell_agents, index_t agent_index,
							   const base_agent_data& base_data, real_t x, real_t y, real_t z, real_t radius_sq,
							   std::vector<index_t>& neighbors)
{
	for (index_t const other_index : cell_agents)
	{
		if (agent_index == other_index)
			continue;

		real_t const ox = base_data.positions[other_index * 3];
		real_t const oy = base_data.positions[other_index * 3 + 1];
		real_t const oz = base_data.positions[other_index * 3 + 2];

		real_t const dist_sq = (x - ox) * (x - ox) + (y - oy) * (y - oy) + (z - oz) * (z - oz);
		if (dist_sq <= radius_sq)
		{
			neighbors.push_back(other_index);
		}
	}
}

} // namespace

uniform_grid_spatial_index::uniform_grid_spatial_index(real_t cell_size) : cell_size(cell_size) {}

void uniform_grid_spatial_index::build(const environment& env)
{
	grid.clear();
	const auto& agents = *env.agents;
	const auto& mech_data_ptr = std::get<std::unique_ptr<agent_data>>(agents.agent_datas);
	const auto& base_data = mech_data_ptr->base_data;
	index_t const count = agents.size();

	for (index_t i = 0; i < count; ++i)
	{
		real_t const x = base_data.positions[i * 3];
		real_t const y = base_data.positions[i * 3 + 1];
		real_t const z = base_data.positions[i * 3 + 2];

		grid_key const key { .x = static_cast<int>(std::floor(x / cell_size)),
							 .y = static_cast<int>(std::floor(y / cell_size)),
							 .z = static_cast<int>(std::floor(z / cell_size)) };
		grid[key].push_back(i);
	}
}

std::vector<index_t> uniform_grid_spatial_index::query_neighbors(const environment& env, index_t agent_index,
																 real_t radius) const
{
	std::vector<index_t> neighbors;
	const auto& agents = *env.agents;
	const auto& mech_data_ptr = std::get<std::unique_ptr<agent_data>>(agents.agent_datas);
	const auto& base_data = mech_data_ptr->base_data;

	real_t const x = base_data.positions[agent_index * 3];
	real_t const y = base_data.positions[agent_index * 3 + 1];
	real_t const z = base_data.positions[agent_index * 3 + 2];

	int const cx = static_cast<int>(std::floor(x / cell_size));
	int const cy = static_cast<int>(std::floor(y / cell_size));
	int const cz = static_cast<int>(std::floor(z / cell_size));

	int const search_radius = static_cast<int>(std::ceil(radius / cell_size));
	real_t const radius_sq = radius * radius;

	for (int dx = -search_radius; dx <= search_radius; ++dx)
	{
		for (int dy = -search_radius; dy <= search_radius; ++dy)
		{
			for (int dz = -search_radius; dz <= search_radius; ++dz)
			{
				grid_key const key { .x = cx + dx, .y = cy + dy, .z = cz + dz };
				auto it = grid.find(key);
				if (it != grid.end())
				{
					collect_neighbors_in_cell(it->second, agent_index, base_data, x, y, z, radius_sq, neighbors);
				}
			}
		}
	}
	return neighbors;
}

} // namespace physicore::mechanics::micromechanics
