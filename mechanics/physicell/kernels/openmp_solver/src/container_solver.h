#pragma once

#include <mutex>

#include <common/generic_agent_solver.h>
#include <physicell/environment.h>

namespace physicore::mechanics::physicell::kernels::openmp_solver {

class container_solver : private generic_agent_solver<mechanical_agent>
{
	std::mutex container_mutex;

public:
	void update_cell_container_for_phenotype(environment& e);
};

} // namespace physicore::mechanics::physicell::kernels::openmp_solver
