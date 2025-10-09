#pragma once

#include "../../../include/solver.h"
#include "bulk_solver.h"
#include "cell_solver.h"
#include "diffusion_solver.h"
#include "dirichlet_solver.h"

namespace physicore::biofvm::kernels::openmp_solver {

class openmp_solver : public solver
{
	bool initialized = false;

	bulk_solver b_solver;
	cell_solver c_solver;
	diffusion_solver d_solver;
	[[no_unique_address]] dirichlet_solver dir_solver;

public:
	void solve(microenvironment& m, index_t iterations) override;
	real_t get_substrate_density(index_t s, index_t x, index_t y, index_t z) const override;
};

} // namespace physicore::biofvm::kernels::openmp_solver
