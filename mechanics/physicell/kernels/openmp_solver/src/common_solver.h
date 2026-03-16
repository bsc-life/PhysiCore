#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

#include <common/cartesian_mesh.h>

namespace physicore::mechanics::physicell::kernels::openmp_solver {

class common_solver
{
public:
	using voxel_pos_t = std::array<index_t, 3>;

	static std::size_t get_mesh_index(const voxel_pos_t& position, const cartesian_mesh& mesh);

	static voxel_pos_t get_mesh_position(const real_t* position, const cartesian_mesh& mesh);

	template <typename func_t>
	static void for_each_in_mech_neighborhood_symmetric(const cartesian_mesh& mesh,
														const std::vector<std::vector<index_t>>& cells_in_voxels,
														const voxel_pos_t& position, index_t i, func_t&& f)
	{
		for_each_in_mech_neighborhood(mesh, cells_in_voxels, position, i, [&](index_t cell_idx) {
			if (cell_idx > i)
				f(cell_idx);
		});
	}

	template <typename func_t>
	static void for_each_in_mech_neighborhood(const cartesian_mesh& mesh,
											  const std::vector<std::vector<index_t>>& cells_in_voxels,
											  const voxel_pos_t& position, index_t i, func_t&& f)
	{
		assert(cells_in_voxels.size() == mesh.voxel_count());

		const sindex_t grid_x = static_cast<sindex_t>(mesh.grid_shape[0]);
		const sindex_t grid_y = static_cast<sindex_t>(mesh.grid_shape[1]);
		const sindex_t grid_z = static_cast<sindex_t>(mesh.grid_shape[2]);

		for (sindex_t z = static_cast<sindex_t>(position[2]) - 1; z <= static_cast<sindex_t>(position[2]) + 1; z++)
		{
			if (z < 0 || z >= grid_z)
				continue;

			for (sindex_t y = static_cast<sindex_t>(position[1]) - 1; y <= static_cast<sindex_t>(position[1]) + 1; y++)
			{
				if (y < 0 || y >= grid_y)
					continue;

				for (sindex_t x = static_cast<sindex_t>(position[0]) - 1; x <= static_cast<sindex_t>(position[0]) + 1;
					 x++)
				{
					if (x < 0 || x >= grid_x)
						continue;

					const auto voxel_idx =
						mesh.linearize(static_cast<index_t>(x), static_cast<index_t>(y), static_cast<index_t>(z));
					for (const auto cell_idx : cells_in_voxels[voxel_idx])
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
