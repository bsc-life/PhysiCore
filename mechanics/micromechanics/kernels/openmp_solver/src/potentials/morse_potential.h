#pragma once

#include <micromechanics/potential_interface.h>
#include <micromechanics/simulation_parameters.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Morse potential for soft cell-cell interactions.
 *
 * The Morse potential provides a smooth repulsion-attraction curve:
 *   V(r) = D * (exp(2*a*(1-r/r0)) - 2*exp(a*(1-r/r0)))
 *
 * Where:
 *   D = potential well depth (derived from stiffness and scaling_factor)
 *   a = scaling_factor
 *   r0 = equilibrium_distance
 *
 * Force is derivative of potential with sign for attraction/repulsion.
 *
 * Recycled from legacy_src/src/solver/host/kelvin_voigt_solver.cpp (solve_pair_inter)
 */
class morse_potential : public potential_interface
{
	interaction_config config_;

public:
	explicit morse_potential(interaction_config config);

	void calculate_pairwise_force(const environment& env, index_t agent_i, index_t agent_j, real_t distance, real_t dx,
								  real_t dy, real_t dz, real_t& force_out) const override;

	std::string name() const override;

	real_t max_interaction_distance(const environment& env, index_t agent_i) const override;
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
