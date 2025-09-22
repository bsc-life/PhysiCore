#include "dirichlet_solver.h"

#include "diffusion_solver.h"
#include "omp_helper.h"

using namespace physicore;
using namespace physicore::biofvm;
using namespace physicore::biofvm::kernels::cpu;

auto fix_dims(const index_t* voxel_index, index_t dims)
{
	if (dims == 1)
		return noarr::fix<'x', 'y', 'z'>(voxel_index[0], 0, 0);
	else if (dims == 2)
		return noarr::fix<'x'>(voxel_index[0]) ^ noarr::fix<'y'>(voxel_index[1]) ^ noarr::fix<'z'>(0);
	else if (dims == 3)
		return noarr::fix<'x'>(voxel_index[0]) ^ noarr::fix<'y'>(voxel_index[1]) ^ noarr::fix<'z'>(voxel_index[2]);
	return noarr::fix<'x', 'y', 'z'>(0, 0, 0);
}

void solve_interior(const auto dens_l, real_t* __restrict__ substrate_densities,
					const index_t* __restrict__ dirichlet_voxels, const real_t* __restrict__ dirichlet_values,
					const bool* __restrict__ dirichlet_conditions, index_t substrates_count,
					index_t dirichlet_voxels_count, index_t dims)
{
	if (dirichlet_voxels_count == 0)
		return;

#pragma omp for
	for (index_t voxel_idx = 0; voxel_idx < dirichlet_voxels_count; ++voxel_idx)
	{
		auto subs_l = dens_l ^ fix_dims(dirichlet_voxels + dims * voxel_idx, dims);

		for (index_t s = 0; s < substrates_count; ++s)
		{
			if (dirichlet_conditions[voxel_idx * substrates_count + s])
				(subs_l | noarr::get_at<'s'>(substrate_densities, s)) =
					dirichlet_values[voxel_idx * substrates_count + s];
		}
	}
}

template <typename density_layout_t>
void solve_boundary(real_t* __restrict__ substrate_densities, const real_t* __restrict__ dirichlet_values,
					const bool* __restrict__ dirichlet_conditions, const density_layout_t dens_l)
{
	if (dirichlet_values == nullptr)
		return;

	omp_trav_for_each_no_parallel(noarr::traverser(dens_l), [=](auto state) {
		auto s = noarr::get_index<'s'>(state);

		if (dirichlet_conditions[s])
			(dens_l | noarr::get_at(substrate_densities, state)) = dirichlet_values[s];
	});
}

void solve_boundaries(const auto dens_l, real_t* __restrict__ substrate_densities, microenvironment& m)
{
	solve_boundary(substrate_densities, m.dirichlet_min_boundary_values[0].get(),
				   m.dirichlet_min_boundary_conditions[0].get(), dens_l ^ noarr::fix<'x'>(noarr::lit<0>));
	solve_boundary(substrate_densities, m.dirichlet_max_boundary_values[0].get(),
				   m.dirichlet_max_boundary_conditions[0].get(), dens_l ^ noarr::fix<'x'>(m.mesh.grid_shape[0] - 1));

	if (m.mesh.dims > 1)
	{
		solve_boundary(substrate_densities, m.dirichlet_min_boundary_values[1].get(),
					   m.dirichlet_min_boundary_conditions[1].get(), dens_l ^ noarr::fix<'y'>(noarr::lit<0>));
		solve_boundary(substrate_densities, m.dirichlet_max_boundary_values[1].get(),
					   m.dirichlet_max_boundary_conditions[1].get(),
					   dens_l ^ noarr::fix<'y'>(m.mesh.grid_shape[1] - 1));
	}

	if (m.mesh.dims > 2)
	{
		solve_boundary(substrate_densities, m.dirichlet_min_boundary_values[2].get(),
					   m.dirichlet_min_boundary_conditions[2].get(), dens_l ^ noarr::fix<'z'>(noarr::lit<0>));
		solve_boundary(substrate_densities, m.dirichlet_max_boundary_values[2].get(),
					   m.dirichlet_max_boundary_conditions[2].get(),
					   dens_l ^ noarr::fix<'z'>(m.mesh.grid_shape[2] - 1));
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
