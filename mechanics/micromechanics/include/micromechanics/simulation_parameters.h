#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>

#include <common/types.h>

namespace physicore::mechanics::micromechanics {

/**
 * @brief Configuration for a specific interaction potential.
 *
 * This struct holds all parameters needed to configure any potential type.
 * Each potential type uses a subset of these parameters.
 */
struct interaction_config
{
	/// Name of the potential to use: "standard", "morse", "kelvin_voigt"
	std::string potential_name = "standard";

	// === Common parameters (used by most potentials) ===

	/// Adhesion strength [force/distance^2]
	real_t adhesion_strength = 0.4;

	/// Repulsion strength [force/distance^2]
	real_t repulsion_strength = 10.0;

	/// Relative distance for maximum adhesion (multiple of cell radius)
	real_t relative_maximum_adhesion_distance = 1.25;

	// === Kelvin-Voigt specific ===

	/// Spring constant for Kelvin-Voigt model
	real_t spring_constant = 1.0;

	/// Damping coefficient for Kelvin-Voigt model
	real_t damping_coefficient = 0.1;

	// === Morse potential specific ===

	/// Scaling factor for Morse potential
	real_t morse_scaling_factor = 1.0;

	/// Equilibrium distance for Morse potential
	real_t morse_equilibrium_distance = 1.0;

	/// Stiffness for Morse potential
	real_t morse_stiffness = 1.0;
};

/**
 * @brief Hash function for agent type pairs.
 *
 * Used for the unordered_map that stores per-type-pair interaction configs.
 */
struct type_pair_hash
{
	std::size_t operator()(const std::pair<std::uint8_t, std::uint8_t>& p) const noexcept
	{
		return std::hash<std::uint16_t> {}((static_cast<std::uint16_t>(p.first) << 8) | p.second);
	}
};

/**
 * @brief Global simulation parameters for the micromechanics system.
 *
 * This struct holds both:
 * - Type-based interaction configurations (which potential + parameters for each type pair)
 * - Global simulation settings (timestep, solver choice, feature flags)
 */
struct simulation_parameters
{
	// === Type-based interaction rules ===

	/// Interaction configuration per agent type pair: (type_a, type_b) -> config
	std::unordered_map<std::pair<std::uint8_t, std::uint8_t>, interaction_config, type_pair_hash> interactions;

	/// Default configuration used when no specific pair config exists
	interaction_config default_interaction;

	// === Global solver settings ===

	/// Name of the solver backend to use (e.g., "openmp_solver")
	std::string solver_name = "openmp_solver";

	/// Mechanics timestep [time units]
	real_t mechanics_timestep = 0.1;

	/// Number of dimensions (2 or 3)
	index_t dims = 3;

	// === Feature flags ===

	bool enable_motility = true;
	bool enable_basement_membrane = false;
	bool enable_spring_attachments = false;

	// === Basement membrane / boundary settings ===

	/// Repulsion strength at domain boundaries
	real_t cell_BM_repulsion_strength = 10.0;

	// === Spring attachment settings ===

	/// Maximum number of spring attachments per cell
	index_t maximum_number_of_attachments = 12;

	/// Spring constant for attachments
	real_t attachment_elastic_constant = 0.01;

	/// Rate at which attachments form [1/time]
	real_t attachment_rate = 0.0;

	/// Rate at which attachments break [1/time]
	real_t detachment_rate = 0.0;

	// === Convenience methods ===

	/**
	 * @brief Add symmetric interaction configuration for a type pair.
	 *
	 * Adds config for both (type_a, type_b) and (type_b, type_a).
	 */
	void add_interaction(std::uint8_t type_a, std::uint8_t type_b, const interaction_config& config)
	{
		interactions[{ type_a, type_b }] = config;
		if (type_a != type_b)
		{
			interactions[{ type_b, type_a }] = config;
		}
	}

	/**
	 * @brief Get interaction config for a type pair.
	 *
	 * Returns specific config if defined, otherwise returns default_interaction.
	 */
	const interaction_config& get_interaction(std::uint8_t type_a, std::uint8_t type_b) const
	{
		auto it = interactions.find({ type_a, type_b });
		if (it != interactions.end())
		{
			return it->second;
		}
		return default_interaction;
	}

	/**
	 * @brief Convenience: set up single-type simulation.
	 *
	 * Sets the (0, 0) interaction to the given config.
	 */
	void set_single_type_interaction(const interaction_config& config) { interactions[{ 0, 0 }] = config; }
};

} // namespace physicore::mechanics::micromechanics
