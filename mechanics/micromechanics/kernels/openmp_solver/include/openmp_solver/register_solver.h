#pragma once

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Attach the OpenMP solver to the global solver registry.
 *
 * This function should be called once at startup to register the
 * OpenMP-based micromechanics solver with the solver_registry.
 */
void attach_to_registry();

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver

