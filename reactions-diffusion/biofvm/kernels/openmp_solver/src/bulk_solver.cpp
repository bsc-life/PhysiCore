#include "bulk_solver.h"

#include <noarr/structures_extended.hpp>

using namespace physicore;
using namespace physicore::biofvm;
using namespace physicore::biofvm::kernels::openmp_solver;

void bulk_solver::initialize(microenvironment& m) { fnc = std::move(m.bulk_fnc); }

template <typename density_layout_t>
void solve_single(real_t* HWY_RESTRICT densities, bulk_functor* fnc, real_t time_step, const density_layout_t dens_l)
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

					const real_t S = fnc->supply_rates(s, x, y, z);
					const real_t U = fnc->uptake_rates(s, x, y, z);
					const real_t T = fnc->supply_target_densities(s, x, y, z);
					real_t& D = dens_l | noarr::get_at(densities, idx);

					D = (D + time_step * S * T) / (1 + time_step * (U + S));
				}
			}
		}
	}
}

void bulk_solver::solve(const microenvironment& m, diffusion_solver& d_solver)
{
	if (fnc)
		solve_single(d_solver.get_substrates_pointer(), fnc.get(), m.diffusion_timestep,
					 d_solver.get_substrates_layout());
}
