#include "thrust_solver/register_solver.h"

#include <biofvm/solver_registry.h>

#include "namespace_config.h"
#include "thrust_solver.h"

void physicore::reactions_diffusion::biofvm::kernels::PHYSICORE_THRUST_SOLVER_NAMESPACE::attach_to_registry()
{
	static const physicore::reactions_diffusion::biofvm::registry_adder<thrust_solver> thrust_solver_adder(
		PHYSICORE_THRUST_SOLVER_REGISTRY_NAME);
}
