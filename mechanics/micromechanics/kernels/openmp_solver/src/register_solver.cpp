#include <micromechanics/solver_registry.h>
#include <openmp_solver/register_solver.h>

#include "openmp_solver.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void attach_to_registry()
{
	// Register the OpenMP solver with the global registry
	solver_registry::instance().register_factory("openmp_solver", []() { return std::make_unique<openmp_solver>(); });
}

// Static registration - solver is registered when the library is loaded
static bool registered = []() {
	attach_to_registry();
	return true;
}();

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
