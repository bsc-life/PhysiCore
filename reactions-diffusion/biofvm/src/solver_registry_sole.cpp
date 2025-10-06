#include <cassert>

#include "solver_registry.h"

using namespace physicore::biofvm;

bool solver_registry::register_factory(std::string solver_name, solver_factory_func_t&& f)
{
	auto [it, emplaced] = factory_registry.try_emplace(std::move(solver_name), std::move(f));

	return emplaced;
}

std::unique_ptr<solver> solver_registry::get(const std::string& solver_name)
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
