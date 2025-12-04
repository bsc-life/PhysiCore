#pragma once

#include <micromechanics/potential_interface.h>
#include <micromechanics/simulation_parameters.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Kelvin-Voigt spring-damper model for cell-cell interactions.
 *
 * Combines spring and damper in parallel:
 *   F = k * (r - r0) + gamma * (v_rel . n)
 *
 * Where:
 *   k = spring constant
 *   r0 = rest length (typically 2 * radius for adjacent cells)
 *   gamma = damping coefficient
 *   v_rel = relative velocity between agents
 *   n = unit normal from agent_i to agent_j
 *
 * Recycled from legacy_src/src/solver/host/kelvin_voigt_solver.cpp (solve_pair_intra)
 *
 * Note: This potential requires velocity information, which is handled specially
 * by the force_solver when using Kelvin-Voigt.
 */
class kelvin_voigt_potential : public potential_interface
{
	interaction_config config_;

public:
	explicit kelvin_voigt_potential(interaction_config  config);

	void calculate_pairwise_force(const environment& env, index_t agent_i, index_t agent_j, real_t distance, real_t dx,
								  real_t dy, real_t dz, real_t& force_out) const override;

	std::string name() const override;

	real_t max_interaction_distance(const environment& env, index_t agent_i) const override;
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
