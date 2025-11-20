#pragma once

#include <thrust/device_new.h>

#include "../../../include/microenvironment.h"
#include "bulk_functor.h"
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

class bulk_solver
{
	thrust::device_ptr<device_bulk_functor> func;

public:
	bulk_solver() = default;
	bulk_solver(const bulk_solver&) = delete;
	bulk_solver& operator=(const bulk_solver&) = delete;

	void initialize(thrust::device_ptr<device_bulk_functor> func);

	void solve(const microenvironment& m, diffusion_solver& d_solver);

	~bulk_solver();
};

} // namespace physicore::biofvm::kernels::thrust_solver
