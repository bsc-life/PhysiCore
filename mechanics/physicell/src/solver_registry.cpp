#include "physicell/solver_registry.h"

#include <cassert>

namespace physicore::mechanics::physicell {

bool solver_registry::register_factory(std::string solver_name, solver_factory_func_t&& f)
{
	auto [it, emplaced] = factory_registry.try_emplace(std::move(solver_name), std::move(f));
	return emplaced;
}

solver_ptr solver_registry::get(const std::string& solver_name)
{
	if (!factory_registry.contains(solver_name))
	{
		assert(false);
		return nullptr;
	}

	return factory_registry[solver_name]();
}

solver_registry& solver_registry::instance()
{
	static solver_registry r;
	return r;
}

} // namespace physicore::mechanics::physicell

