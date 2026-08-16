#pragma once

#include <physicell/solver.h>

namespace physicore::mechanics::physicell::kernels::openmp_solver {

class openmp_solver : public solver
{
	bool initialized = false;

public:
	void initialize(environment& e) override;
	void solve(environment& e, index_t iterations) override;
};

} // namespace physicore::mechanics::physicell::kernels::openmp_solver
