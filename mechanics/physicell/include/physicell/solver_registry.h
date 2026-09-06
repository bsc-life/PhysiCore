#pragma once

#include <common/registry.h>

#include "solver.h"

namespace physicore::mechanics::physicell {

class solver_registry : public physicore::factory_registry<solver>
{
public:
	static solver_registry& instance();
};

template <typename T>
concept derived_from_solver = std::derived_from<T, solver>;

template <derived_from_solver SolverT>
using registry_adder = physicore::generic_registry_adder<SolverT, solver_registry>;

} // namespace physicore::mechanics::physicell
