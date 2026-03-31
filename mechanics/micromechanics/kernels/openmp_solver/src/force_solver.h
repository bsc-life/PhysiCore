#pragma once

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include <common/types.h>
#include <micromechanics/potential_interface.h>
#include <micromechanics/simulation_parameters.h>

namespace physicore::mechanics::micromechanics {
class environment;
}

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Handles force calculations between agents using configurable potentials.
 *
 * The force solver implements type-based interactions: different agent type pairs
 * can use different potentials (Standard, Morse, Kelvin-Voigt). It creates and
 * caches potential instances based on the simulation configuration.
 */
class force_solver
{
	bool initialized_ = false;

	/// Cached potentials for each agent type pair: (type_a, type_b) -> potential
	std::unordered_map<std::pair<std::uint8_t, std::uint8_t>, std::unique_ptr<potential_interface>, type_pair_hash>
		interaction_potentials_;

	/// Default potential for unconfigured type pairs
	std::unique_ptr<potential_interface> default_potential_;

public:
	void initialize(environment& e);

	/**
	 * @brief Calculate forces for all agents.
	 *
	 * Iterates over all agent pairs within interaction distance and
	 * applies the appropriate potential based on their types.
	 * Forces are accumulated in the agent_data.forces array.
	 */
	void calculate_forces(environment& e);

private:
	/**
	 * @brief Create a potential instance from configuration.
	 */
	static std::unique_ptr<potential_interface> create_potential(const interaction_config& config);

	/**
	 * @brief Get the potential for a type pair.
	 */
	const potential_interface& get_potential(std::uint8_t type_a, std::uint8_t type_b) const;

	/**
	 * @brief Clear all force accumulators.
	 */
	static void clear_forces(environment& e);
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
