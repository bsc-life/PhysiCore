#pragma once

#include <micromechanics/potential_interface.h>
#include <micromechanics/simulation_parameters.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Morse potential for soft cell-cell interactions (r²/r₀² form).
 *
 * Uses a modified Morse potential with squared-distance exponent:
 *   V(r) = D * [ exp(2·a·(1 − r²/r₀²)) − 2·exp(a·(1 − r²/r₀²)) ]
 *
 * Where:
 *   D  = potential well depth = (k · r₀²) / (8 · a²)
 *   a  = scaling_factor  (controls well width)
 *   r₀ = equilibrium_distance
 *   k  = stiffness
 *
 * The force (−dV/dr) evaluates to:
 *   F(r) = (4·a·r·D / r₀²) · [ exp(2P) − exp(P) ]
 *   with  P = a · (1 − r²/r₀²)
 *
 * The r²/r₀² form is intentional: it gives a steeper repulsive wall
 * and smoother adhesion tail than the standard r/r₀ Morse potential.
 */
class morse_potential : public potential_interface
{
	interaction_config config_;

public:
	explicit morse_potential(interaction_config config);

	real_t calculate_pairwise_force(const environment& env, index_t agent_i, index_t agent_j, real_t distance,
									real_t dx, real_t dy, real_t dz) const override;

	std::string name() const override;

	real_t max_interaction_distance(const environment& env, index_t agent_i) const override;
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
