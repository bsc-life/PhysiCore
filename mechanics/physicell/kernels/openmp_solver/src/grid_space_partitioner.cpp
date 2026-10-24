#include "grid_space_partitioner.h"

#include <noarr/structures/extra/shortcuts.hpp>

namespace physicore::mechanics::physicell::kernels::openmp_solver {

void grid_space_partitioner::initialize(index_t voxel_size, const cartesian_mesh& microenv_mesh)
{
	partitioning_mesh = cartesian_mesh(microenv_mesh.dims, microenv_mesh.bounding_box_mins,
									   microenv_mesh.bounding_box_maxs, { voxel_size, voxel_size, voxel_size });
	const index_t voxels_count = partitioning_mesh.voxel_count();
	agents_in_voxels_sizes = std::make_unique<std::atomic<index_t>[]>(voxels_count);
	agents_in_voxels = std::make_unique<std::vector<index_t>[]>(voxels_count);
}

template <>
index_t grid_space_partitioner::get_mesh_index<1>(std::array<index_t, 3> point) const
{
	auto mesh_l = noarr::scalar<uint8_t>() ^ noarr::vectors<'x'>(partitioning_mesh.grid_shape[0]);
	return mesh_l | noarr::offset<'x'>(point[0]);
}

template <>
index_t grid_space_partitioner::get_mesh_index<2>(std::array<index_t, 3> point) const
{
	auto mesh_l = noarr::scalar<uint8_t>()
				  ^ noarr::vectors<'x', 'y'>(partitioning_mesh.grid_shape[0], partitioning_mesh.grid_shape[1]);
	return mesh_l | noarr::offset<'x', 'y'>(point[0], point[1]);
}

template <>
index_t grid_space_partitioner::get_mesh_index<3>(std::array<index_t, 3> point) const
{
	auto mesh_l = noarr::scalar<uint8_t>()
				  ^ noarr::vectors<'x', 'y', 'z'>(partitioning_mesh.grid_shape[0], partitioning_mesh.grid_shape[1],
												  partitioning_mesh.grid_shape[2]);
	return mesh_l | noarr::offset<'x', 'y', 'z'>(point[0], point[1], point[2]);
}

index_t grid_space_partitioner::get_mesh_index(const real_t* position) const
{
	if (partitioning_mesh.dims == 1)
	{
		auto voxel_pos = partitioning_mesh.voxel_position({ position, 1 });

		return get_mesh_index<1>(voxel_pos);
	}
	if (partitioning_mesh.dims == 2)
	{
		auto voxel_pos = partitioning_mesh.voxel_position({ position, 2 });

		return get_mesh_index<2>(voxel_pos);
	}
	{
		auto voxel_pos = partitioning_mesh.voxel_position({ position, 3 });

		return get_mesh_index<3>(voxel_pos);
	}
}

void grid_space_partitioner::update_partitioning(const real_t* PHYSICORE_RESTRICT positions, index_t agents_count)
{
#pragma omp for
	for (std::size_t i = 0; i < partitioning_mesh.voxel_count(); i++)
	{
		agents_in_voxels[i].clear();
		agents_in_voxels_sizes[i].store(0, std::memory_order_relaxed);
	}

	// first we count how many cells are in each voxel
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		agents_in_voxels_sizes[get_mesh_index(positions + i * partitioning_mesh.dims)].fetch_add(
			1, std::memory_order_relaxed);
	}

	// second we allocate memory for each voxel
#pragma omp for
	for (std::size_t i = 0; i < partitioning_mesh.voxel_count(); i++)
	{
		agents_in_voxels[i].resize(agents_in_voxels_sizes[i].load(std::memory_order_relaxed));
	}

	// third we assign cells to voxels
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		auto mech_idx = get_mesh_index(positions + i * partitioning_mesh.dims);

		auto in_voxel_index = agents_in_voxels_sizes[mech_idx].fetch_sub(1, std::memory_order_relaxed) - 1;

		agents_in_voxels[mech_idx][in_voxel_index] = i;
	}
}

} // namespace physicore::mechanics::physicell::kernels::openmp_solver
