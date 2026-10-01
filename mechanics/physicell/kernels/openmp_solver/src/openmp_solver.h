#pragma once

#include <physicell/solver.h>

#include "position_solver.h"

namespace physicore::mechanics::physicell::kernels::openmp_solver {

class openmp_solver : public solver
{
	bool initialized = false;

	position_solver mechanics_position_solver;

public:
	void initialize(environment& e) override;
	void solve(environment& e, index_t iterations) override;
	migration_bias_func_ptr create_migration_bias_functor(environment& e, migration_bias_type type) override;
};

} // namespace physicore::mechanics::physicell::kernels::openmp_solver
