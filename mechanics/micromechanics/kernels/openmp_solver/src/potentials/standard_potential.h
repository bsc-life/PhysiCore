#pragma once

#include <micromechanics/potential_interface.h>
#include <micromechanics/simulation_parameters.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Standard PhysiCell-style cell-cell interaction potential.
 *
 * Implements the classic PhysiCell repulsion-adhesion model:
 * - Repulsion: quadratic when distance < sum_of_radii
 * - Adhesion: quadratic when distance < relative_max_adhesion_distance * sum_of_radii
 *
 * Force = repulsion_strength * (1 - d/R)^2 - adhesion_strength * (1 - d/R_adh)^2
 *
 * Recycled from legacy_src/src/solver/host/standard_position_solver.cpp
 */
class standard_potential : public potential_interface
{
	interaction_config config_;

public:
	explicit standard_potential(const interaction_config& config);

	void calculate_pairwise_force(const environment& env, index_t agent_i, index_t agent_j, real_t distance, real_t dx,
								  real_t dy, real_t dz, real_t& force_out) const override;

	std::string name() const override;

	real_t max_interaction_distance(const environment& env, index_t agent_i) const override;
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
