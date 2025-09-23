#include "cpu_solver/register_solver.h"

#include "../../../include/solver_registry.h"
#include "cpu_solver.h"

void physicore::biofvm::kernels::cpu::attach_to_registry()
{
	static physicore::biofvm::registry_adder<cpu_solver> cpu_solver_adder("cpu_solver");
}
