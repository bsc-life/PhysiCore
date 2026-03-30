#pragma once

#include <common/types.h>

namespace physicore::mechanics::micromechanics {
class environment;
}

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Handles position integration.
 *
 * Integrates agent positions using the Adams-Bashforth 2nd order method,
 * which provides improved stability compared to simple Euler integration.
 */
class position_solver
{
	bool initialized_ = false;

public:
	void initialize(environment& e);

	/**
	 * @brief Update agent positions using accumulated forces/velocities.
	 *
	 * Uses Adams-Bashforth 2nd order integration:
	 *   x_new = x_old + dt * (1.5 * v_new - 0.5 * v_old)
	 *
	 * Also updates previous_velocities for next timestep.
	 */
	void update_positions(environment& e);
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
