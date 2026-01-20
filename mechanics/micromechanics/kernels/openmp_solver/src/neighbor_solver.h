#pragma once

#include <common/types.h>

namespace physicore::mechanics::micromechanics {
class environment;
}

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Handles neighbor finding and spatial indexing.
 *
 * Responsible for building and querying spatial data structures
 * to efficiently find nearby agents for force calculations.
 */
class neighbor_solver
{
	bool initialized_ = false;

public:
	void initialize(environment& e);

	/**
	 * @brief Update neighbor lists for all agents.
	 *
	 * Rebuilds the spatial index and populates neighbor lists
	 * for each agent based on interaction distances.
	 */
	static void update_neighbors(environment& e);
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
