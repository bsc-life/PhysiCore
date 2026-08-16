#pragma once

#include <array>
#include <cassert>
#include <span>

#include "types.h"

namespace physicore {

struct cartesian_mesh
{
	index_t dims; // 1 or 2 or 3

	std::array<sindex_t, 3> bounding_box_mins; // [x_min, y_min, z_min]
	std::array<sindex_t, 3> bounding_box_maxs; // [x_max, y_max, z_max]

	std::array<index_t, 3> voxel_shape; // [dx, dy, dz]
	std::array<index_t, 3> grid_shape;	// [x_size, y_size, z_size]

	constexpr cartesian_mesh(index_t dims, std::array<sindex_t, 3> bounding_box_mins,
							 std::array<sindex_t, 3> bounding_box_maxs, std::array<index_t, 3> voxel_shape) noexcept
		: dims(dims),
		  bounding_box_mins(bounding_box_mins),
		  bounding_box_maxs(bounding_box_maxs),
		  voxel_shape(voxel_shape),
		  grid_shape({ 1, 1, 1 })
	{
		if (dims >= 1)
		{
			grid_shape[0] = (bounding_box_maxs[0] - bounding_box_mins[0] + voxel_shape[0] - 1) / voxel_shape[0];
		}
		if (dims >= 2)
		{
			grid_shape[1] = (bounding_box_maxs[1] - bounding_box_mins[1] + voxel_shape[1] - 1) / voxel_shape[1];
		}
		if (dims >= 3)
		{
			grid_shape[2] = (bounding_box_maxs[2] - bounding_box_mins[2] + voxel_shape[2] - 1) / voxel_shape[2];
		}
	}


	constexpr std::size_t voxel_count() const noexcept
	{
		return (std::size_t)grid_shape[0] * (std::size_t)grid_shape[1] * (std::size_t)grid_shape[2];
	}

	constexpr index_t voxel_volume() const noexcept { return voxel_shape[0] * voxel_shape[1] * voxel_shape[2]; }

	constexpr std::size_t linearize(index_t x, index_t y, index_t z) const noexcept
	{
		return x + y * grid_shape[0] + z * grid_shape[0] * grid_shape[1];
	}

	constexpr std::array<index_t, 3> voxel_position(std::span<const real_t> position) const noexcept
	{
		assert(position.size() == (size_t)dims);
		assert(position[0] <= bounding_box_maxs[0] && position[0] >= bounding_box_mins[0]);
		if (dims >= 2)
			assert(position[1] <= bounding_box_maxs[1] && position[1] >= bounding_box_mins[1]);
		if (dims >= 3)
			assert(position[2] <= bounding_box_maxs[2] && position[2] >= bounding_box_mins[2]);

		switch (position.size())
		{
			case 1:
				return { (index_t)((position[0] - (real_t)bounding_box_mins[0]) / (real_t)voxel_shape[0]), 0, 0 };
			case 2:
				return { (index_t)((position[0] - (real_t)bounding_box_mins[0]) / (real_t)voxel_shape[0]),
						 (index_t)((position[1] - (real_t)bounding_box_mins[1]) / (real_t)voxel_shape[1]), 0 };
			case 3:
				return { (index_t)((position[0] - (real_t)bounding_box_mins[0]) / (real_t)voxel_shape[0]),
						 (index_t)((position[1] - (real_t)bounding_box_mins[1]) / (real_t)voxel_shape[1]),
						 (index_t)((position[2] - (real_t)bounding_box_mins[2]) / (real_t)voxel_shape[2]) };
			default:
				assert(false); // Should never reach here
				return { 0, 0, 0 };
		}
	}

	constexpr std::array<real_t, 3> voxel_center(std::array<index_t, 3> position) const noexcept
	{
		assert(position[0] < grid_shape[0]);
		assert(position[1] < grid_shape[1]);
		assert(position[2] < grid_shape[2]);

		return { (real_t)(position[0] * voxel_shape[0] + bounding_box_mins[0]) + ((real_t)voxel_shape[0] / (real_t)2.0),
				 (real_t)(position[1] * voxel_shape[1] + bounding_box_mins[1]) + ((real_t)voxel_shape[1] / (real_t)2.0),
				 (real_t)(position[2] * voxel_shape[2] + bounding_box_mins[2])
					 + ((real_t)voxel_shape[2] / (real_t)2.0) };
	}
};

} // namespace physicore
