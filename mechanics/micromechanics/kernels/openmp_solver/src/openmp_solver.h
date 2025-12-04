#pragma once

#include <micromechanics/solver.h>

#include "basement_membrane_solver.h"
#include "force_solver.h"
#include "motility_solver.h"
#include "neighbor_solver.h"
#include "position_solver.h"
#include "spring_solver.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief OpenMP-parallelized micromechanics solver.
 *
 * This solver implements all micromechanics operations using OpenMP
 * for shared-memory parallelization. It delegates to specialized
 * sub-solvers for each operation type.
 *
 * The solver is registered with the solver_registry as "openmp_solver"
 * and can be instantiated via:
 *   solver_registry::instance().get("openmp_solver")
 */
class openmp_solver : public solver
{
	bool initialized_ = false;

	neighbor_solver n_solver_;
	force_solver f_solver_;
	motility_solver m_solver_;
	basement_membrane_solver bm_solver_;
	spring_solver s_solver_;
	position_solver p_solver_;

public:
	void initialize(environment& e) override;
	void update_cell_neighbors(environment& e) override;
	void update_cell_forces(environment& e) override;
	void calculate_cell_data(environment& e) override;
	void update_motility(environment& e) override;
	void update_basement_membrane_interactions(environment& e) override;
	void update_spring_attachments(environment& e) override;
	void update_positions(environment& e) override;
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
