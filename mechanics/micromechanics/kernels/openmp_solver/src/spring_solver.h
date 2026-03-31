#pragma once

#include <common/types.h>

namespace physicore::mechanics::micromechanics {
class environment;
}

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Handles spring attachment dynamics.
 *
 * Manages the formation, breakage, and force calculation of
 * spring attachments between cells. Springs can form between
 * nearby cells and break based on detachment rates.
 */
class spring_solver
{
	bool initialized_ = false;

public:
	void initialize(environment& e);

	/**
	 * @brief Update spring attachments and calculate spring forces.
	 *
	 * 1. Mark springs for detachment based on detachment_rate
	 * 2. Form new springs between nearby cells based on attachment_rate
	 * 3. Calculate spring contraction forces
	 */
	static void update_spring_attachments(environment& e);
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
