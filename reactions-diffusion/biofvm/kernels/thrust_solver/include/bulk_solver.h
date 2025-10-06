#pragma once

#include <thrust/device_new.h>

#include "../../../include/microenvironment.h"
#include "diffusion_solver.h"

/*
Performs bulk supply and uptake of substrates.

For each bulk (voxel), the following is performed:
S = supply_rate_f(m, voxel_idx)
U = uptake_rate_f(m, voxel_idx)
T = supply_target_densities_f(m, voxel_idx)
D = (D + dt*S*T)/(1 + dt*(U+S))
where D is a voxel substrate density vector
*/

namespace physicore::biofvm::kernels::thrust_solver {

struct bulk_functor
{
	PHYSICORE_THRUST_DEVICE_FN virtual real_t supply_rates(index_t s, index_t x, index_t y, index_t z) = 0;
	PHYSICORE_THRUST_DEVICE_FN virtual real_t uptake_rates(index_t s, index_t x, index_t y, index_t z) = 0;
	PHYSICORE_THRUST_DEVICE_FN virtual real_t supply_target_densities(index_t s, index_t x, index_t y, index_t z) = 0;
	PHYSICORE_THRUST_DEVICE_FN virtual ~bulk_functor() {}
};

class bulk_solver
{
	thrust::device_ptr<bulk_functor> func;


public:
	void initialize(const microenvironment& m);

	template <typename FuncType>
	void initialize()
	{
		func = thrust::device_new<FuncType>();
	}

	void solve(const microenvironment& m, diffusion_solver& d_solver);

	~bulk_solver();
};

} // namespace physicore::biofvm::kernels::thrust_solver
