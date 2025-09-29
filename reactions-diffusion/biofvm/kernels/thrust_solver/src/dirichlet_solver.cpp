#include "dirichlet_solver.h"

#include <hwy/highway.h>

using namespace physicore;
using namespace physicore::biofvm;
using namespace physicore::biofvm::kernels::thrust_solver;

static auto fix_dims(const index_t* voxel_index, index_t dims)
{
	if (dims == 1)
		return noarr::fix<'x', 'y', 'z'>(voxel_index[0], 0, 0);
	else if (dims == 2)
		return noarr::fix<'x'>(voxel_index[0]) ^ noarr::fix<'y'>(voxel_index[1]) ^ noarr::fix<'z'>(0);
	else if (dims == 3)
		return noarr::fix<'x'>(voxel_index[0]) ^ noarr::fix<'y'>(voxel_index[1]) ^ noarr::fix<'z'>(voxel_index[2]);
	return noarr::fix<'x', 'y', 'z'>(0, 0, 0);
}

static void solve_interior(const auto dens_l, real_t* HWY_RESTRICT substrate_densities,
						   const index_t* HWY_RESTRICT dirichlet_voxels, const real_t* HWY_RESTRICT dirichlet_values,
						   const bool* HWY_RESTRICT dirichlet_conditions, index_t substrates_count,
						   index_t dirichlet_voxels_count, index_t dims)
{
	if (dirichlet_voxels_count == 0)
		return;

	thrust::for_each(thrust::make_counting_iterator<index_t>(0), thrust::make_counting_iterator(dirichlet_voxels_count),
					 [=](index_t voxel_idx) {
						 auto subs_l = dens_l ^ fix_dims(dirichlet_voxels + dims * voxel_idx, dims);

						 for (index_t s = 0; s < substrates_count; ++s)
						 {
							 if (dirichlet_conditions[voxel_idx * substrates_count + s])
								 (subs_l | noarr::get_at<'s'>(substrate_densities, s)) =
									 dirichlet_values[voxel_idx * substrates_count + s];
						 }
					 });
}

static void solve_boundaries(const auto dens_l, real_t* HWY_RESTRICT substrate_densities, microenvironment& m)
{
	index_t s_len = m.substrates_count;
	index_t x_len = m.mesh.grid_shape[0];
	index_t y_len = m.mesh.grid_shape[1];
	index_t z_len = m.mesh.grid_shape[2];

	std::size_t dirichlet_x_side_work = (std::size_t)s_len * y_len * z_len;

	// x min
	if (m.dirichlet_min_boundary_values[0])
		thrust::for_each(thrust::make_counting_iterator<std::size_t>(0),
						 thrust::make_counting_iterator(dirichlet_x_side_work), [=, &m](std::size_t voxel_idx) {
							 const index_t s = voxel_idx % s_len;
							 voxel_idx /= s_len;
							 const index_t y = voxel_idx % y_len;
							 const index_t z = voxel_idx / y_len;

							 if (m.dirichlet_min_boundary_conditions[0][s])
								 (dens_l | noarr::get_at<'s', 'x', 'y', 'z'>(substrate_densities, s, 0, y, z)) =
									 m.dirichlet_min_boundary_values[0][s];
						 });

	// x max
	if (m.dirichlet_max_boundary_values[0])
		thrust::for_each(thrust::make_counting_iterator<std::size_t>(0),
						 thrust::make_counting_iterator(dirichlet_x_side_work), [=, &m](std::size_t voxel_idx) {
							 const index_t s = voxel_idx % s_len;
							 voxel_idx /= s_len;
							 const index_t y = voxel_idx % y_len;
							 const index_t z = voxel_idx / y_len;

							 if (m.dirichlet_max_boundary_conditions[0][s])
								 (dens_l | noarr::get_at<'s', 'x', 'y', 'z'>(substrate_densities, s, x_len - 1, y, z)) =
									 m.dirichlet_max_boundary_values[0][s];
						 });


	if (m.mesh.dims > 1)
	{
		std::size_t dirichlet_y_side_work = (std::size_t)s_len * x_len * z_len;

		// y min
		if (m.dirichlet_min_boundary_values[1])
			thrust::for_each(thrust::make_counting_iterator<std::size_t>(0),
							 thrust::make_counting_iterator(dirichlet_y_side_work), [=, &m](std::size_t voxel_idx) {
								 const index_t s = voxel_idx % s_len;
								 voxel_idx /= s_len;
								 const index_t x = voxel_idx % x_len;
								 const index_t z = voxel_idx / x_len;

								 if (m.dirichlet_min_boundary_conditions[1][s])
									 (dens_l | noarr::get_at<'s', 'x', 'y', 'z'>(substrate_densities, s, x, 0, z)) =
										 m.dirichlet_min_boundary_values[1][s];
							 });

		// y max
		if (m.dirichlet_max_boundary_values[1])
			thrust::for_each(thrust::make_counting_iterator<std::size_t>(0),
							 thrust::make_counting_iterator(dirichlet_y_side_work), [=, &m](std::size_t voxel_idx) {
								 const index_t s = voxel_idx % s_len;
								 voxel_idx /= s_len;
								 const index_t x = voxel_idx % x_len;
								 const index_t z = voxel_idx / x_len;

								 if (m.dirichlet_max_boundary_conditions[1][s])
									 (dens_l
									  | noarr::get_at<'s', 'x', 'y', 'z'>(substrate_densities, s, x, y_len - 1, z)) =
										 m.dirichlet_max_boundary_values[1][s];
							 });
	}

	if (m.mesh.dims > 2)
	{
		std::size_t dirichlet_z_side_work = (std::size_t)s_len * x_len * y_len;

		// z min
		if (m.dirichlet_min_boundary_values[2])
			thrust::for_each(thrust::make_counting_iterator<std::size_t>(0),
							 thrust::make_counting_iterator(dirichlet_z_side_work), [=, &m](std::size_t voxel_idx) {
								 const index_t s = voxel_idx % s_len;
								 voxel_idx /= s_len;
								 const index_t x = voxel_idx % x_len;
								 const index_t y = voxel_idx / x_len;

								 if (m.dirichlet_min_boundary_conditions[2][s])
									 (dens_l | noarr::get_at<'s', 'x', 'y', 'z'>(substrate_densities, s, x, y, 0)) =
										 m.dirichlet_min_boundary_values[2][s];
							 });

		// z max
		if (m.dirichlet_max_boundary_values[2])
			thrust::for_each(thrust::make_counting_iterator<std::size_t>(0),
							 thrust::make_counting_iterator(dirichlet_z_side_work), [=, &m](std::size_t voxel_idx) {
								 const index_t s = voxel_idx % s_len;
								 voxel_idx /= s_len;
								 const index_t x = voxel_idx % x_len;
								 const index_t y = voxel_idx / x_len;

								 if (m.dirichlet_max_boundary_conditions[2][s])
									 (dens_l
									  | noarr::get_at<'s', 'x', 'y', 'z'>(substrate_densities, s, x, y, z_len - 1)) =
										 m.dirichlet_max_boundary_values[2][s];
							 });
	}
}

void dirichlet_solver::solve(microenvironment& m, diffusion_solver& d_solver)
{
	solve_boundaries(d_solver.get_substrates_layout(), d_solver.get_substrates_pointer(), m);
	solve_interior(d_solver.get_substrates_layout(), d_solver.get_substrates_pointer(),
				   m.dirichlet_interior_voxels.get(), m.dirichlet_interior_values.get(),
				   m.dirichlet_interior_conditions.get(), m.substrates_count, m.dirichlet_interior_voxels_count,
				   m.mesh.dims);
}
