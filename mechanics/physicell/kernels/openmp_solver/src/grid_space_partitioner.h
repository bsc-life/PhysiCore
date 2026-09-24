#pragma once

#include <atomic>
#include <memory>
#include <vector>

#include <common/mesh.h>

namespace physicore::mechanics::physicell::kernels::openmp_solver {

class grid_space_partitioner
{
	cartesian_mesh partitioning_mesh {};

	std::unique_ptr<std::atomic<index_t>[]> agents_in_voxels_sizes;
	std::unique_ptr<std::vector<index_t>[]> agents_in_voxels;

	template <index_t dims>
	index_t get_mesh_index(std::array<index_t, 3> point) const;
	index_t get_mesh_index(const real_t* position) const;

public:
	void initialize(index_t voxel_size, const cartesian_mesh& microenv_mesh);

	void update_partitioning(const real_t* positions, index_t agents_count);

	template <index_t dims, typename func_t>
	void for_each_in_neighborhood(const real_t* agent_position, index_t i, func_t f) const
	{
		auto position = partitioning_mesh.voxel_position({ agent_position, dims });
		for (sindex_t z = -1; z <= 1; z++)
		{
			if (position[2] + z >= partitioning_mesh.grid_shape[2] || (sindex_t)position[2] + z < 0)
				continue;

			for (sindex_t y = -1; y <= 1; y++)
			{
				if (position[1] + y >= partitioning_mesh.grid_shape[1] || (sindex_t)position[1] + y < 0)
					continue;

				for (sindex_t x = -1; x <= 1; x++)
				{
					if (position[0] + x >= partitioning_mesh.grid_shape[0] || (sindex_t)position[0] + x < 0)
						continue;

					index_t voxel_index {};

					if constexpr (dims == 1)
					{
						voxel_index = get_mesh_index<1>({ position[0] + x, 0, 0 });
					}
					else if constexpr (dims == 2)
					{
						voxel_index = get_mesh_index<2>({ position[0] + x, position[1] + y, 0 });
					}
					else
					{
						voxel_index = get_mesh_index<3>({ position[0] + x, position[1] + y, position[2] + z });
					}

					for (auto& cell_idx : agents_in_voxels[voxel_index])
					{
						if (i != cell_idx)
							f(cell_idx);
					}
				}
			}
		}
	}
};

} // namespace physicore::mechanics::physicell::kernels::openmp_solver
