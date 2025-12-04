#pragma once

#include <concepts>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "solver.h"

namespace physicore::mechanics::micromechanics {

/**
 * @brief Registry for micromechanics solver backends.
 *
 * Following the BioFVM pattern, this registry allows runtime selection of
 * solver implementations. Backends self-register via the registry_adder template
 * at static initialization time.
 *
 * Usage:
 *   auto solver = solver_registry::instance().get("openmp_solver");
 */
class solver_registry
{
public:
	using solver_factory_func_t = std::function<solver_ptr()>;
	using registry_map_t = std::unordered_map<std::string, solver_factory_func_t>;

	registry_map_t factory_registry;

	/**
	 * @brief Register a solver factory function.
	 * @param solver_name Unique name for the solver
	 * @param f Factory function that creates solver instances
	 * @return true if registration succeeded, false if name already exists
	 */
	bool register_factory(std::string solver_name, solver_factory_func_t&& f);

	/**
	 * @brief Create a solver instance by name.
	 * @param solver_name Name of the solver to create
	 * @return Unique pointer to solver, or nullptr if not found
	 */
	solver_ptr get(const std::string& solver_name);

	/**
	 * @brief Check if a solver is available.
	 * @param solver_name Name of the solver to check
	 * @return true if solver is registered
	 */
	bool is_available(const std::string& solver_name) const;

	/**
	 * @brief Get list of all registered solver names.
	 * @return Vector of solver names
	 */
	std::vector<std::string> available_solvers() const;

	/**
	 * @brief Get the global solver registry instance.
	 * @return Reference to singleton registry
	 */
	static solver_registry& instance();
};

/// Concept requiring a type to derive from solver
template <typename T>
concept derived_from_solver = std::derived_from<T, solver>;

/**
 * @brief Helper class for automatic solver registration.
 *
 * Create a static instance of this class to register a solver at startup.
 * Example:
 *   static registry_adder<my_solver> my_solver_adder("my_solver");
 */
template <derived_from_solver SolverT>
struct registry_adder
{
	explicit registry_adder(std::string solver_name)
	{
		solver_registry::instance().register_factory(std::move(solver_name),
													 []() { return std::make_unique<SolverT>(); });
	}
};

} // namespace physicore::mechanics::micromechanics
