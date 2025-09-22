#include "bulk_solver.h"

#include <noarr/structures_extended.hpp>

using namespace physicore;
using namespace physicore::biofvm::kernels::cpu;

void bulk_solver::initialize(microenvironment& m)
{
	supply_rate_f_ = m.supply_rate_func;
	uptake_rate_f_ = m.uptake_rate_func;
	supply_target_densities_f_ = m.supply_target_densities_func;

	supply_rates_ = std::make_unique<real_t[]>(m.substrates_count);
	uptake_rates_ = std::make_unique<real_t[]>(m.substrates_count);
	supply_target_densities_ = std::make_unique<real_t[]>(m.substrates_count);
}

template <typename density_layout_t>
void solve_single(real_t* HWY_RESTRICT densities, const real_t* HWY_RESTRICT supply_rates,
				  const real_t* HWY_RESTRICT uptake_rates, const real_t* HWY_RESTRICT supply_target_densities,
				  real_t time_step, const density_layout_t dens_l)
{
	const index_t s_dim = dens_l | noarr::get_length<'s'>();

	for (index_t s = 0; s < s_dim; s++)
	{
		(dens_l | noarr::get_at<'s'>(densities, s)) =
			((dens_l | noarr::get_at<'s'>(densities, s)) + time_step * supply_rates[s] * supply_target_densities[s])
			/ (1 + time_step * (uptake_rates[s] + supply_rates[s]));
	}
}

void bulk_solver::solve(microenvironment& m, diffusion_solver& d_solver)
{
#pragma omp for collapse(3)
	for (index_t z = 0; z < m.mesh.grid_shape[2]; z++)
	{
		for (index_t y = 0; y < m.mesh.grid_shape[1]; y++)
		{
			for (index_t x = 0; x < m.mesh.grid_shape[0]; x++)
			{
				supply_rate_f_(m, { x, y, z }, supply_rates_.get());
				uptake_rate_f_(m, { x, y, z }, uptake_rates_.get());
				supply_target_densities_f_(m, { x, y, z }, supply_target_densities_.get());

				solve_single(d_solver.get_substrates_pointer(), supply_rates_.get(), uptake_rates_.get(),
							 supply_target_densities_.get(), m.diffusion_timestep,
							 d_solver.get_substrates_layout() ^ noarr::fix<'x'>(x) ^ noarr::fix<'y'>(y)
								 ^ noarr::fix<'z'>(z));
			}
		}
	}
}
