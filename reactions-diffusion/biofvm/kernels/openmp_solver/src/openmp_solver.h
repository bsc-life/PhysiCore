#pragma once

#include <biofvm/solver.h>

#include "bulk_solver.h"
#include "cell_solver.h"
#include "diffusion_solver.h"

namespace physicore::biofvm::kernels::openmp_solver {

class openmp_solver : public solver
{
	bool initialized = false;
	bool recompute_cells = true;

	bulk_solver b_solver;
	cell_solver c_solver;
	diffusion_solver d_solver;

public:
	void initialize(microenvironment& m) override;
	void solve(microenvironment& m, index_t iterations) override;
	real_t get_substrate_density(index_t s, index_t x, index_t y, index_t z) const override;
	real_t& get_substrate_density(index_t s, index_t x, index_t y, index_t z) override;
	void reinitialize_dirichlet(microenvironment& m) override;
	void recompute_positional_data(microenvironment& m) override;
};

} // namespace physicore::biofvm::kernels::openmp_solver
