#include "spring_solver.h"

#include <cmath>
#include <random>
#include <tuple>

#include <micromechanics/agent_container.h>
#include <micromechanics/agent_data.h>
#include <micromechanics/environment.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void spring_solver::initialize(environment& /*e*/)
{
	if (initialized_)
		return;

	initialized_ = true;
}

void spring_solver::update_spring_attachments(environment& e)
{
	if (!e.params.enable_spring_attachments)
		return;

	auto& agents = *e.agents;
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);
	auto& base_data = mech_data.base_data;
	index_t const count = agents.size();
	real_t const dt = e.timestep;

	real_t const elastic_constant = e.params.attachment_elastic_constant;
	real_t const attachment_rate = e.params.attachment_rate;
	real_t const detachment_rate = e.params.detachment_rate;

	// Note: Spring attachment management is complex and requires careful
	// synchronization. For now, we implement a simplified version that
	// calculates forces for existing springs.
	//
	// Full implementation would need:
	// 1. Thread-safe spring creation/deletion
	// 2. Neighbor queries for potential new attachments
	// 3. Proper handling of maximum_number_of_attachments

#pragma omp parallel for
	for (index_t i = 0; i < count; ++i)
	{
		if (!mech_data.is_movable[i])
			continue;

		// Get spring attachments for this agent
		const auto& attachments = mech_data.spring_attachments[i];

		for (index_t const j : attachments)
		{
			if (j >= count)
				continue; // Invalid attachment

			// Calculate spring force
			real_t const dx = base_data.positions[j * 3] - base_data.positions[i * 3];
			real_t const dy = base_data.positions[j * 3 + 1] - base_data.positions[i * 3 + 1];
			real_t const dz = base_data.positions[j * 3 + 2] - base_data.positions[i * 3 + 2];
			real_t const distance = std::sqrt(dx * dx + dy * dy + dz * dz);

			if (distance < 1e-16)
				continue;

			// Rest length is sum of radii (touching)
			real_t const rest_length = mech_data.radii[i] + mech_data.radii[j];

			// Spring force: F = k * (distance - rest_length)
			// Positive when stretched (pulls together), negative when compressed
			real_t const force_mag = elastic_constant * (distance - rest_length);

			// Normalize direction
			real_t const nx = dx / distance;
			real_t const ny = dy / distance;
			real_t const nz = dz / distance;

			// Apply force (pulls i toward j when stretched)
			mech_data.forces[i * 3] += force_mag * nx;
			mech_data.forces[i * 3 + 1] += force_mag * ny;
			mech_data.forces[i * 3 + 2] += force_mag * nz;
		}
	}

	// TODO: Implement spring formation and breakage
	// This would require:
	// 1. Random number generation for attachment/detachment probabilities
	// 2. Thread-safe modification of spring_attachments vectors
	// 3. Neighbor queries to find potential attachment partners
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
