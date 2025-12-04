#include <cassert>

#include <micromechanics/solver_registry.h>

namespace physicore::mechanics::micromechanics {

bool solver_registry::register_factory(std::string solver_name, solver_factory_func_t&& f)
{
	auto [it, emplaced] = factory_registry.try_emplace(std::move(solver_name), std::move(f));
	return emplaced;
}

solver_ptr solver_registry::get(const std::string& solver_name)
{
	if (!factory_registry.contains(solver_name))
	{
		assert(false && "Unknown solver requested");
		return nullptr;
	}

	return factory_registry[solver_name]();
}

bool solver_registry::is_available(const std::string& solver_name) const
{
	return factory_registry.contains(solver_name);
}

std::vector<std::string> solver_registry::available_solvers() const
{
	std::vector<std::string> names;
	names.reserve(factory_registry.size());
	for (const auto& [name, factory] : factory_registry)
	{
		names.push_back(name);
	}
	return names;
}

solver_registry& solver_registry::instance()
{
	static solver_registry r;
	return r;
}

} // namespace physicore::mechanics::micromechanics
