#pragma once

#include <biofvm/microenvironment.h>

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

namespace physicore::biofvm::kernels::openmp_solver {

class bulk_solver
{
	std::unique_ptr<bulk_functor> fnc;

public:
	void initialize(microenvironment& m);

	void solve(const microenvironment& m, diffusion_solver& d_solver);
};

} // namespace physicore::biofvm::kernels::openmp_solver
