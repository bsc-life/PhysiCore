#pragma once

#include <common/types.h>

namespace physicore::mechanics::micromechanics {
class environment;
}

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Handles cell motility forces.
 *
 * Implements random walk with persistence, migration bias,
 * and optional chemotaxis for cell motility.
 */
class motility_solver
{
	bool initialized_ = false;

public:
	void initialize(environment& e);

	/**
	 * @brief Update motility forces for all motile agents.
	 *
	 * For each motile agent:
	 * - With probability dt/persistence_time, update motility direction
	 * - Apply random walk with optional bias direction
	 * - Scale by migration speed
	 * - Add to velocity accumulator
	 */
	static void update_motility(environment& e);
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
