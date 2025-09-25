#include "bulk_solver.h"

#include <noarr/structures_extended.hpp>

using namespace physicore;
using namespace physicore::biofvm::kernels::cpu;

void bulk_solver::initialize(const microenvironment& m)
{
	supply_rate_f_ = m.supply_rate_func;
	uptake_rate_f_ = m.uptake_rate_func;
	supply_target_densities_f_ = m.supply_target_densities_func;
}

template <typename density_layout_t>
void solve_single(real_t* HWY_RESTRICT densities, auto&& supply_rates, auto&& uptake_rates,
				  auto&& supply_target_densities, real_t time_step, const density_layout_t dens_l)
{
	const index_t x_dim = dens_l | noarr::get_length<'x'>();
	const index_t y_dim = dens_l | noarr::get_length<'y'>();
	const index_t z_dim = dens_l | noarr::get_length<'z'>();
	const index_t s_dim = dens_l | noarr::get_length<'s'>();

#pragma omp for collapse(3)
	for (index_t z = 0; z < z_dim; z++)
	{
		for (index_t y = 0; y < y_dim; y++)
		{
			for (index_t x = 0; x < x_dim; x++)
			{
				for (index_t s = 0; s < s_dim; s++)
				{
					auto idx = noarr::idx<'s', 'x', 'y', 'z'>(s, x, y, z);

					const real_t S = supply_rates(s, x, y, z);
					const real_t U = uptake_rates(s, x, y, z);
					const real_t T = supply_target_densities(s, x, y, z);
					real_t& D = dens_l | noarr::get_at(densities, idx);

					D = (D + time_step * S * T) / (1 + time_step * (U + S));
				}
			}
		}
	}
}

void bulk_solver::solve(const microenvironment& m, diffusion_solver& d_solver)
{
	solve_single(d_solver.get_substrates_pointer(), supply_rate_f_, uptake_rate_f_, supply_target_densities_f_,
				 m.diffusion_timestep, d_solver.get_substrates_layout());
}
