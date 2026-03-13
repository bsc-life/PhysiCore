#pragma once

#include <common/registry.h>

#include "solver.h"

namespace physicore::mechanics::physicell {

class solver_registry : public physicore::common::factory_registry<solver>
{
public:
	using base_t = physicore::common::factory_registry<solver>;
	using solver_factory_func_t = typename base_t::factory_func_t;
	using registry_map_t = typename base_t::map_t;

	static solver_registry& instance();
};

template <typename SolverT>
using registry_adder = physicore::common::registry_adder<SolverT, solver_registry, solver>;

} // namespace physicore::mechanics::physicell
