#pragma once

#include <array>
#include <span>

#include "types.h"

namespace physicore {

/**
 * @brief Uniform Cartesian mesh for spatial domain discretization.
 *
 * Provides axis-aligned structured grid with uniform spacing.
 * Used by both reaction-diffusion and mechanics modules for domain partitioning.
 *
 * Supports 1D, 2D, and 3D domains with configurable voxel sizes and bounding boxes.
 * Enables efficient spatial queries: position-to-voxel mapping and voxel linearization.
 *
 * @see reactions-diffusion/biofvm for diffusion solver usage
 * @see mechanics/physicell for mechanics engine usage
 */
class cartesian_mesh
{
public:
	/// @brief Constructor
	/// @param dims Number of spatial dimensions (1, 2, or 3)
	/// @param bounding_box_mins Minimum coordinates of domain bounds
	/// @param bounding_box_maxs Maximum coordinates of domain bounds
	/// @param voxel_shape Size of each voxel in each dimension
	explicit cartesian_mesh(index_t dims, std::array<sindex_t, 3> bounding_box_mins,
							std::array<sindex_t, 3> bounding_box_maxs, std::array<index_t, 3> voxel_shape);

	/// @brief Get total number of voxels in the mesh
	[[nodiscard]] std::size_t voxel_count() const;

	/// @brief Get volume of a single voxel
	[[nodiscard]] index_t voxel_volume() const;

	/// @brief Find voxel indices containing the given position
	/// @param position Spatial coordinates (must match dims)
	/// @return Voxel indices in x, y, z order
	[[nodiscard]] std::array<index_t, 3> voxel_position(std::span<const real_t> position) const;

	/// @brief Get the center position of a voxel
	/// @param position Voxel indices
	/// @return Center coordinates in x, y, z order
	[[nodiscard]] std::array<real_t, 3> voxel_center(std::array<index_t, 3> position) const;

	/// @brief Convert 3D voxel indices to linear index
	/// @param x, y, z Voxel indices
	/// @return Linear index for use in flat arrays
	[[nodiscard]] std::size_t linearize(index_t x, index_t y, index_t z) const;

	// Member data (public for direct access, following biofvm pattern)
	index_t dims;							   ///< Number of spatial dimensions
	std::array<sindex_t, 3> bounding_box_mins; ///< Minimum domain coordinates
	std::array<sindex_t, 3> bounding_box_maxs; ///< Maximum domain coordinates
	std::array<index_t, 3> voxel_shape;		   ///< Size of each voxel
	std::array<index_t, 3> grid_shape;		   ///< Number of voxels per dimension
};

} // namespace physicore
