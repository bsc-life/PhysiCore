#include "force_solver.h"

#include <algorithm>
#include <cmath>
#include <tuple>

#include <micromechanics/agent_container.h>
#include <micromechanics/agent_data.h>
#include <micromechanics/environment.h>
#include <micromechanics/spatial_index.h>

#include "potentials/kelvin_voigt_potential.h"
#include "potentials/morse_potential.h"
#include "potentials/standard_potential.h"

#ifdef _OPENMP
	#include <omp.h>
#endif

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void force_solver::initialize(environment& e)
{
	if (initialized_)
		return;

	// Create potentials for all configured type pairs
	for (const auto& [type_pair, config] : e.params.interactions)
	{
		interaction_potentials_[type_pair] = create_potential(config);
	}

	// Create default potential
	default_potential_ = create_potential(e.params.default_interaction);

	initialized_ = true;
}

std::unique_ptr<potential_interface> force_solver::create_potential(const interaction_config& config)
{
	if (config.potential_name == "morse")
	{
		return std::make_unique<morse_potential>(config);
	}
	if (config.potential_name == "kelvin_voigt")
	{
		return std::make_unique<kelvin_voigt_potential>(config);
	}
	else // default to standard
	{
		return std::make_unique<standard_potential>(config);
	}
}

const potential_interface& force_solver::get_potential(std::uint8_t type_a, std::uint8_t type_b) const
{
	auto it = interaction_potentials_.find({ type_a, type_b });
	if (it != interaction_potentials_.end())
	{
		return *it->second;
	}
	return *default_potential_;
}

void force_solver::clear_forces(environment& e)
{
	auto& agents = *e.agents;
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);

	std::ranges::fill(mech_data.forces, 0.0);
}

void force_solver::calculate_forces(environment& e)
{
	auto& agents = *e.agents;
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);
	auto& base_data = mech_data.base_data;
	index_t const count = agents.size();

	// Clear force accumulators
	clear_forces(e);

#pragma omp parallel for
	for (index_t i = 0; i < count; ++i)
	{
		if (!mech_data.is_movable[i])
			continue;

		// Get agent type from agent_data
		std::uint8_t const type_i = mech_data.agent_types[i];

		// Determine max interaction distance from default potential
		real_t const max_dist = default_potential_->max_interaction_distance(e, i);

		// Query neighbors within interaction distance
		auto neighbors = e.index->query_neighbors(e, i, max_dist);

		for (index_t const j : neighbors)
		{
			std::uint8_t const type_j = mech_data.agent_types[j];

			// Calculate position difference (j - i, so positive = j is ahead of i)
			// Legacy uses this convention: position_difference points from i to j
			real_t const dx = base_data.positions[j * 3] - base_data.positions[i * 3];
			real_t const dy = base_data.positions[j * 3 + 1] - base_data.positions[i * 3 + 1];
			real_t const dz = base_data.positions[j * 3 + 2] - base_data.positions[i * 3 + 2];
			real_t distance = std::sqrt(dx * dx + dy * dy + dz * dz);

			// Avoid division by zero (legacy uses 0.00001)
			distance = std::max(distance, 0.00001);

			// Get appropriate potential for this type pair
			const auto& potential = get_potential(type_i, type_j);

			// Calculate force coefficient (already divided by distance in potential)
			// Legacy: force = (repulsion - adhesion) / distance
			// Then: velocity += position_difference * force
			real_t force_coeff = 0.0;
			potential.calculate_pairwise_force(e, i, j, distance, dx, dy, dz, force_coeff);

			// Accumulate force using position difference
			// dx points from i toward j (dx = pos_j - pos_i)
			// force_coeff > 0 means net repulsion: i should move AWAY from j (opposite to dx)
			// force_coeff < 0 means net adhesion: i should move TOWARD j (same as dx)
			// So we use NEGATIVE force_coeff to get correct direction:
			//   repulsion (force_coeff > 0) → -force_coeff * dx → i moves away from j
			//   adhesion (force_coeff < 0) → -force_coeff * dx → i moves toward j
			mech_data.forces[i * 3] -= force_coeff * dx;
			mech_data.forces[i * 3 + 1] -= force_coeff * dy;
			mech_data.forces[i * 3 + 2] -= force_coeff * dz;
		}
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
