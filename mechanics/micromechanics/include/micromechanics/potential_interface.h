#pragma once

#include <string>

#include <common/types.h>

namespace physicore::mechanics::micromechanics {

class environment;
struct interaction_config;

/**
 * @brief Abstract interface for pairwise interaction potentials.
 *
 * Potentials define the force calculation between pairs of agents.
 * Different potential types (Standard, Morse, Kelvin-Voigt) implement
 * different force models while sharing a common interface.
 *
 * Potentials are configured per agent-type pair, allowing heterogeneous
 * interactions in the simulation (e.g., cell type A-A uses Morse,
 * A-B uses Standard).
 */
class potential_interface
{
public:
	/**
	 * @brief Calculate force between two agents.
	 *
	 * @param env The simulation environment
	 * @param agent_i Index of first agent
	 * @param agent_j Index of second agent
	 * @param distance Pre-calculated distance between agents
	 * @param dx X-component of normalized direction vector (j - i)
	 * @param dy Y-component of normalized direction vector
	 * @param dz Z-component of normalized direction vector
	 * @param force_out Output: force magnitude (positive = repulsion, negative = attraction)
	 */
	virtual void calculate_pairwise_force(const environment& env, index_t agent_i, index_t agent_j, real_t distance,
										  real_t dx, real_t dy, real_t dz, real_t& force_out) const = 0;

	/// Get the name of this potential type
	virtual std::string name() const = 0;

	/// Get the maximum interaction distance for this potential
	virtual real_t max_interaction_distance(const environment& env, index_t agent_i) const = 0;

	virtual ~potential_interface() = default;
};

} // namespace physicore::mechanics::micromechanics

