#pragma once

#include <array>
#include <span>

#include <common/types.h>

namespace physicore::biofvm {

struct cartesian_mesh
{
	index_t dims; // 1 or 2 or 3

	std::array<sindex_t, 3> bounding_box_mins; // [x_min, y_min, z_min]
	std::array<sindex_t, 3> bounding_box_maxs; // [x_max, y_max, z_max]

	std::array<index_t, 3> voxel_shape; // [dx, dy, dz]
	std::array<index_t, 3> grid_shape;	// [x_size, y_size, z_size]

	cartesian_mesh(index_t dims, std::array<sindex_t, 3> bounding_box_mins, std::array<sindex_t, 3> bounding_box_maxs,
				   std::array<index_t, 3> voxel_shape);

	std::size_t voxel_count() const;
	index_t voxel_volume() const;

	std::size_t linearize(index_t x, index_t y, index_t z) const;

	std::array<index_t, 3> voxel_position(std::span<const real_t> position) const;

	std::array<real_t, 3> voxel_center(std::array<index_t, 3> position) const;
};

} // namespace physicore::biofvm
