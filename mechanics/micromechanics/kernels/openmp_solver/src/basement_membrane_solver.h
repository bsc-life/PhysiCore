#pragma once

#include <common/types.h>

namespace physicore::mechanics::micromechanics {
class environment;
}

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Handles basement membrane and boundary interactions.
 *
 * Applies repulsive forces when agents approach domain boundaries,
 * preventing them from leaving the simulation domain.
 */
class basement_membrane_solver
{
	bool initialized_ = false;

public:
	void initialize(environment& e);

	/**
	 * @brief Update basement membrane interaction forces.
	 *
	 * For each agent near a domain boundary, applies a repulsive
	 * force to keep the agent inside the domain.
	 */
	void update_interactions(environment& e);
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
