#pragma once

#include <concepts>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

#include "solver.h"

namespace physicore::mechanics::physicell {

class solver_registry
{
public:
	using solver_factory_func_t = std::function<solver_ptr()>;
	using registry_map_t = std::unordered_map<std::string, solver_factory_func_t>;

	registry_map_t factory_registry;

	bool register_factory(std::string solver_name, solver_factory_func_t&& f);

	solver_ptr get(const std::string& solver_name);

	static solver_registry& instance();
};

template <typename T>
concept derived_from_solver = std::derived_from<T, solver>;

template <derived_from_solver SolverT>
struct registry_adder
{
	explicit registry_adder(std::string solver_name)
	{
		solver_registry::instance().register_factory(std::move(solver_name),
													 []() { return std::make_unique<SolverT>(); });
	}
};

} // namespace physicore::mechanics::physicell
