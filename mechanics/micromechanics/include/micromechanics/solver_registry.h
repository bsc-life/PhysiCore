#pragma once

#include <concepts>
#include <string>

#include <common/factory_registry.h>

#include "solver.h"

namespace physicore::mechanics::micromechanics {

/// Registry for micromechanics solver backends.
///
/// This is implemented via the generic `physicore::factory_registry` (in `common/`).
using solver_registry = physicore::factory_registry<solver>;

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
using registry_adder = physicore::factory_registry_adder<solver, SolverT>;

} // namespace physicore::mechanics::micromechanics
