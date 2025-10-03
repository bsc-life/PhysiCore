#include "dirichlet_solver.h"

#include <thrust/execution_policy.h>
#include <thrust/for_each.h>

using namespace physicore;
using namespace physicore::biofvm;
using namespace physicore::biofvm::kernels::thrust_solver;

static constexpr auto fix_dims(const index_t* voxel_index, index_t dims)
{
	if (dims == 1)
		return noarr::fix<'x', 'y', 'z'>(voxel_index[0], 0, 0);
	else if (dims == 2)
		return noarr::fix<'x'>(voxel_index[0]) ^ noarr::fix<'y'>(voxel_index[1]) ^ noarr::fix<'z'>(0);
	else if (dims == 3)
		return noarr::fix<'x'>(voxel_index[0]) ^ noarr::fix<'y'>(voxel_index[1]) ^ noarr::fix<'z'>(voxel_index[2]);
	return noarr::fix<'x', 'y', 'z'>(0, 0, 0);
}

template <typename density_layout_t>
static void solve_interior(const density_layout_t dens_l, real_t* _CCCL_RESTRICT substrate_densities,
						   const index_t* _CCCL_RESTRICT dirichlet_voxels,
						   const real_t* _CCCL_RESTRICT dirichlet_values,
						   const bool* _CCCL_RESTRICT dirichlet_conditions, index_t substrates_count,
						   index_t dirichlet_voxels_count, index_t dims)
{
	if (dirichlet_voxels_count == 0)
		return;

	thrust::for_each(thrust::device, thrust::make_counting_iterator<index_t>(0),
					 thrust::make_counting_iterator(dirichlet_voxels_count),
					 [dens_l, substrate_densities, dirichlet_voxels, dirichlet_values, dirichlet_conditions,
					  substrates_count, dims] PHYSICORE_THRUST_DEVICE_FN(index_t voxel_idx) {
						 auto subs_l = dens_l ^ fix_dims(dirichlet_voxels + dims * voxel_idx, dims);

						 for (index_t s = 0; s < substrates_count; ++s)
						 {
							 if (dirichlet_conditions[voxel_idx * substrates_count + s])
								 (subs_l | noarr::get_at<'s'>(substrate_densities, s)) =
									 dirichlet_values[voxel_idx * substrates_count + s];
						 }
					 });
}

template <typename density_layout_t>
static void solve_boundaries(const density_layout_t dens_l, real_t* _CCCL_RESTRICT substrate_densities,
							 microenvironment& m,
							 std::array<thrust::device_vector<real_t>, 3>& dirichlet_min_boundary_values,
							 std::array<thrust::device_vector<real_t>, 3>& dirichlet_max_boundary_values,
							 std::array<thrust::device_vector<bool>, 3>& dirichlet_min_boundary_conditions,
							 std::array<thrust::device_vector<bool>, 3>& dirichlet_max_boundary_conditions)
{
	index_t s_len = m.substrates_count;
	index_t x_len = m.mesh.grid_shape[0];
	index_t y_len = m.mesh.grid_shape[1];
	index_t z_len = m.mesh.grid_shape[2];

	std::size_t dirichlet_x_side_work = (std::size_t)s_len * y_len * z_len;

	// x min
	if (m.dirichlet_min_boundary_values[0])
		thrust::for_each(
			thrust::device, thrust::make_counting_iterator<std::size_t>(0),
			thrust::make_counting_iterator(dirichlet_x_side_work),
			[dens_l, s_len, y_len, substrate_densities,
			 dirichlet_conditions = dirichlet_min_boundary_conditions[0].data().get(),
			 dirichlet_values =
				 dirichlet_min_boundary_values[0].data().get()] PHYSICORE_THRUST_DEVICE_FN(std::size_t voxel_idx) {
				const index_t s = voxel_idx % s_len;
				voxel_idx /= s_len;
				const index_t y = voxel_idx % y_len;
				const index_t z = voxel_idx / y_len;

				if (dirichlet_conditions[s])
					(dens_l | noarr::get_at<'s', 'x', 'y', 'z'>(substrate_densities, s, 0, y, z)) = dirichlet_values[s];
			});

	// x max
	if (m.dirichlet_max_boundary_values[0])
		thrust::for_each(
			thrust::device, thrust::make_counting_iterator<std::size_t>(0),
			thrust::make_counting_iterator(dirichlet_x_side_work),
			[dens_l, s_len, x_len, y_len, substrate_densities,
			 dirichlet_conditions = dirichlet_max_boundary_conditions[0].data().get(),
			 dirichlet_values =
				 dirichlet_max_boundary_values[0].data().get()] PHYSICORE_THRUST_DEVICE_FN(std::size_t voxel_idx) {
				const index_t s = voxel_idx % s_len;
				voxel_idx /= s_len;
				const index_t y = voxel_idx % y_len;
				const index_t z = voxel_idx / y_len;

				if (dirichlet_conditions[s])
					(dens_l | noarr::get_at<'s', 'x', 'y', 'z'>(substrate_densities, s, x_len - 1, y, z)) =
						dirichlet_values[s];
			});


	if (m.mesh.dims > 1)
	{
		std::size_t dirichlet_y_side_work = (std::size_t)s_len * x_len * z_len;

		// y min
		if (m.dirichlet_min_boundary_values[1])
			thrust::for_each(
				thrust::device, thrust::make_counting_iterator<std::size_t>(0),
				thrust::make_counting_iterator(dirichlet_y_side_work),
				[dens_l, s_len, x_len, substrate_densities,
				 dirichlet_conditions = dirichlet_min_boundary_conditions[1].data().get(),
				 dirichlet_values =
					 dirichlet_min_boundary_values[1].data().get()] PHYSICORE_THRUST_DEVICE_FN(std::size_t voxel_idx) {
					const index_t s = voxel_idx % s_len;
					voxel_idx /= s_len;
					const index_t x = voxel_idx % x_len;
					const index_t z = voxel_idx / x_len;

					if (dirichlet_conditions[s])
						(dens_l | noarr::get_at<'s', 'x', 'y', 'z'>(substrate_densities, s, x, 0, z)) =
							dirichlet_values[s];
				});

		// y max
		if (m.dirichlet_max_boundary_values[1])
			thrust::for_each(
				thrust::device, thrust::make_counting_iterator<std::size_t>(0),
				thrust::make_counting_iterator(dirichlet_y_side_work),
				[dens_l, s_len, x_len, y_len, substrate_densities,
				 dirichlet_conditions = dirichlet_max_boundary_conditions[1].data().get(),
				 dirichlet_values =
					 dirichlet_max_boundary_values[1].data().get()] PHYSICORE_THRUST_DEVICE_FN(std::size_t voxel_idx) {
					const index_t s = voxel_idx % s_len;
					voxel_idx /= s_len;
					const index_t x = voxel_idx % x_len;
					const index_t z = voxel_idx / x_len;

					if (dirichlet_conditions[s])
						(dens_l | noarr::get_at<'s', 'x', 'y', 'z'>(substrate_densities, s, x, y_len - 1, z)) =
							dirichlet_values[s];
				});
	}

	if (m.mesh.dims > 2)
	{
		std::size_t dirichlet_z_side_work = (std::size_t)s_len * x_len * y_len;

		// z min
		if (m.dirichlet_min_boundary_values[2])
			thrust::for_each(
				thrust::device, thrust::make_counting_iterator<std::size_t>(0),
				thrust::make_counting_iterator(dirichlet_z_side_work),
				[dens_l, s_len, x_len, substrate_densities,
				 dirichlet_conditions = dirichlet_min_boundary_conditions[2].data().get(),
				 dirichlet_values =
					 dirichlet_min_boundary_values[2].data().get()] PHYSICORE_THRUST_DEVICE_FN(std::size_t voxel_idx) {
					const index_t s = voxel_idx % s_len;
					voxel_idx /= s_len;
					const index_t x = voxel_idx % x_len;
					const index_t y = voxel_idx / x_len;

					if (dirichlet_conditions[s])
						(dens_l | noarr::get_at<'s', 'x', 'y', 'z'>(substrate_densities, s, x, y, 0)) =
							dirichlet_values[s];
				});

		// z max
		if (m.dirichlet_max_boundary_values[2])
			thrust::for_each(
				thrust::device, thrust::make_counting_iterator<std::size_t>(0),
				thrust::make_counting_iterator(dirichlet_z_side_work),
				[dens_l, s_len, x_len, z_len, substrate_densities,
				 dirichlet_conditions = dirichlet_max_boundary_conditions[2].data().get(),
				 dirichlet_values =
					 dirichlet_max_boundary_values[2].data().get()] PHYSICORE_THRUST_DEVICE_FN(std::size_t voxel_idx) {
					const index_t s = voxel_idx % s_len;
					voxel_idx /= s_len;
					const index_t x = voxel_idx % x_len;
					const index_t y = voxel_idx / x_len;

					if (dirichlet_conditions[s])
						(dens_l | noarr::get_at<'s', 'x', 'y', 'z'>(substrate_densities, s, x, y, z_len - 1)) =
							dirichlet_values[s];
				});
	}
}

void dirichlet_solver::initialize(microenvironment& m)
{
	index_t s_len = m.substrates_count;

	auto copy = [](auto&& host_pointer, auto&& device_vector, index_t len) {
		if (host_pointer == nullptr)
			return;

		device_vector.resize(len);
		device_vector.shrink_to_fit();

		thrust::copy(host_pointer, host_pointer + len, device_vector.begin());
	};

	for (index_t i = 0; i < m.mesh.dims; i++)
	{
		copy(m.dirichlet_min_boundary_values[i].get(), dirichlet_min_boundary_values[i], s_len);
		copy(m.dirichlet_max_boundary_values[i].get(), dirichlet_max_boundary_values[i], s_len);
		copy(m.dirichlet_min_boundary_conditions[i].get(), dirichlet_min_boundary_conditions[i], s_len);
		copy(m.dirichlet_max_boundary_conditions[i].get(), dirichlet_max_boundary_conditions[i], s_len);
	}

	copy(m.dirichlet_interior_voxels.get(), dirichlet_interior_voxels, m.dirichlet_interior_voxels_count * m.mesh.dims);
	copy(m.dirichlet_interior_values.get(), dirichlet_interior_values, m.dirichlet_interior_voxels_count * s_len);
	copy(m.dirichlet_interior_conditions.get(), dirichlet_interior_conditions,
		 m.dirichlet_interior_voxels_count * s_len);
}

void dirichlet_solver::solve(microenvironment& m, diffusion_solver& d_solver)
{
	solve_boundaries(d_solver.get_substrates_layout(), d_solver.get_substrates_pointer().get(), m,
					 dirichlet_min_boundary_values, dirichlet_max_boundary_values, dirichlet_min_boundary_conditions,
					 dirichlet_max_boundary_conditions);

	if (m.dirichlet_interior_voxels_count != 0)
		solve_interior(d_solver.get_substrates_layout(), d_solver.get_substrates_pointer().get(),
					   dirichlet_interior_voxels.data().get(), dirichlet_interior_values.data().get(),
					   dirichlet_interior_conditions.data().get(), m.substrates_count,
					   m.dirichlet_interior_voxels_count, m.mesh.dims);
}
