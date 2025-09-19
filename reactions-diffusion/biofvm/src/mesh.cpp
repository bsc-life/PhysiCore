#include "mesh.h"

#include <cassert>

using namespace physicore;
using namespace physicore::biofvm;

cartesian_mesh::cartesian_mesh(index_t dims, std::array<index_t, 3> bounding_box_mins,
							   std::array<index_t, 3> bounding_box_maxs, std::array<index_t, 3> voxel_shape)
	: dims(dims), bounding_box_mins(bounding_box_mins), bounding_box_maxs(bounding_box_maxs), voxel_shape(voxel_shape)
{
	grid_shape = { 1, 1, 1 };

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

std::size_t cartesian_mesh::voxel_count() const
{
	return (std::size_t)grid_shape[0] * (std::size_t)grid_shape[1] * (std::size_t)grid_shape[2];
}

index_t cartesian_mesh::voxel_volume() const { return voxel_shape[0] * voxel_shape[1] * voxel_shape[2]; }

std::array<index_t, 3> cartesian_mesh::voxel_position(std::span<real_t> position) const
{
	assert(position.size() == (size_t)dims);

	switch (dims)
	{
		case 1:
			return { (index_t)((position[0] - bounding_box_mins[0]) / voxel_shape[0]),
					 (index_t)((position[1] - bounding_box_mins[1]) / voxel_shape[1]) };
		case 2:
			return { (index_t)((position[0] - bounding_box_mins[0]) / voxel_shape[0]),
					 (index_t)((position[1] - bounding_box_mins[1]) / voxel_shape[1]) };
		case 3:
			return { (index_t)((position[0] - bounding_box_mins[0]) / voxel_shape[0]),
					 (index_t)((position[1] - bounding_box_mins[1]) / voxel_shape[1]),
					 (index_t)((position[2] - bounding_box_mins[2]) / voxel_shape[2]) };
	}

	assert(false); // Should never reach here
	return { 0, 0, 0 };
}

std::array<real_t, 3> cartesian_mesh::voxel_center(std::array<index_t, 3> position) const
{
	return { (real_t)(position[0] * voxel_shape[0] + voxel_shape[0] / 2.0 + bounding_box_mins[0]),
			 (real_t)(position[1] * voxel_shape[1] + voxel_shape[1] / 2.0 + bounding_box_mins[1]),
			 (real_t)(position[2] * voxel_shape[2] + voxel_shape[2] / 2.0 + bounding_box_mins[2]) };
}
