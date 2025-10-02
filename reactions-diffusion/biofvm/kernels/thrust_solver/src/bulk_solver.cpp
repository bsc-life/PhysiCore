#include "bulk_solver.h"

#include <cuda/std/functional>
#include <noarr/structures_extended.hpp>
#include <thrust/execution_policy.h>
#include <thrust/for_each.h>

#include "thrust/iterator/counting_iterator.h"

using namespace physicore;
using namespace physicore::biofvm::kernels::thrust_solver;

void bulk_solver::initialize(const microenvironment& m)
{
	supply_rate_f_ = m.supply_rate_func;
	uptake_rate_f_ = m.uptake_rate_func;
	supply_target_densities_f_ = m.supply_target_densities_func;
}

template <typename density_layout_t, typename func_t>
void solve_single(real_t* _CCCL_RESTRICT densities, func_t&& supply_rates, func_t&& uptake_rates,
				  func_t&& supply_target_densities, real_t time_step, const density_layout_t dens_l)
{
	std::size_t n = (std::size_t)(dens_l | noarr::get_length<'s'>()) * (dens_l | noarr::get_length<'x'>())
					* (dens_l | noarr::get_length<'y'>()) * (dens_l | noarr::get_length<'z'>());

	thrust::for_each(
		thrust::host, thrust::make_counting_iterator<std::size_t>(0), thrust::make_counting_iterator<std::size_t>(n),
		[dens_l, densities, supply_rates, uptake_rates, supply_target_densities, time_step](std::size_t i) {
			const index_t x_dim = dens_l | noarr::get_length<'x'>();
			const index_t y_dim = dens_l | noarr::get_length<'y'>();
			const index_t s_dim = dens_l | noarr::get_length<'s'>();

			const index_t s = i % s_dim;
			i /= s_dim;
			const index_t x = i % x_dim;
			i /= x_dim;
			const index_t y = i % y_dim;
			const index_t z = i / y_dim;

			auto idx = noarr::idx<'s', 'x', 'y', 'z'>(s, x, y, z);

			const real_t S = supply_rates(s, x, y, z);
			const real_t U = uptake_rates(s, x, y, z);
			const real_t T = supply_target_densities(s, x, y, z);

			real_t& D = dens_l | noarr::get_at(densities, idx);

			D = (D + time_step * S * T) / (1 + time_step * (U + S));
		});
}

void bulk_solver::solve(const microenvironment& m, diffusion_solver& d_solver)
{
	solve_single(d_solver.get_substrates_pointer(), supply_rate_f_, uptake_rate_f_, supply_target_densities_f_,
				 m.diffusion_timestep, d_solver.get_substrates_layout());
}
