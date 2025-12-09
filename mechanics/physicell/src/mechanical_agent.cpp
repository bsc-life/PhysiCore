#include "physicell/mechanical_agent.h"

#include <algorithm>
#include <cassert>

namespace physicore::mechanics::physicell {

std::unique_ptr<mechanical_agent> mechanical_agent::add(index_t new_id, index_t cell_type,
														mechanical_parameters& parameters, bool is_2D)
{
	return add(new_id, cell_type, parameters, data, is_2D);
}

std::unique_ptr<mechanical_agent> mechanical_agent::add(index_t new_id, index_t cell_type,
														mechanical_parameters& parameters, agent_data& data, bool is_2D)
{
	assert(new_id == data.base_data.agents_count);

	data.base_data.add();
	data.add();

	assert(new_id < data.agents_count);

	auto agent = std::make_unique<mechanical_agent>(new_id, data);
	agent->index = new_id;
	agent->type_index = cell_type;

	const index_t dim_offset = new_id * data.dims;
	std::fill_n(&data.velocities[dim_offset], data.dims, 0.0);
	std::fill_n(&data.previous_velocities[dim_offset], data.dims, 0.0);
	std::fill_n(&data.migration_bias_direction[dim_offset], data.dims, 0.0);
	std::fill_n(&data.motility_vector[dim_offset], data.dims, 0.0);
	std::fill_n(&data.orientation[dim_offset], data.dims, 0.0);

	data.radius[new_id] = 0.0;

	data.cell_definition_indices[new_id] = cell_type;

	// Directly assign scalar values from mechanical_parameters
	data.cell_cell_adhesion_strength[new_id] = parameters.cell_cell_adhesion_strength;
	data.cell_BM_adhesion_strength[new_id] = parameters.cell_BM_adhesion_strength;
	data.cell_cell_repulsion_strength[new_id] = parameters.cell_cell_repulsion_strength;
	data.cell_BM_repulsion_strength[new_id] = parameters.cell_BM_repulsion_strength;

	// Copy cell_adhesion_affinities vector directly
	const index_t type_count = data.agent_types_count;
	for (index_t other = 0; other < type_count && other < parameters.cell_adhesion_affinity.size(); ++other)
	{
		data.cell_adhesion_affinities[new_id * type_count + other] = parameters.cell_adhesion_affinity[other];
	}

	data.relative_maximum_adhesion_distance[new_id] = parameters.relative_maximum_adhesion_distance;
	data.maximum_number_of_attachments[new_id] = parameters.maximum_number_of_attachments;
	data.attachment_elastic_constant[new_id] = parameters.attachment_elastic_coefficient;
	data.attachment_rate[new_id] = parameters.attachment_rate;
	data.detachment_rate[new_id] = parameters.detachment_rate;

	data.is_motile[new_id] = static_cast<std::uint8_t>(parameters.is_movable);
	data.is_movable[new_id] = static_cast<std::uint8_t>(parameters.is_movable);
	data.persistence_time[new_id] = parameters.motility_persistence_time;
	data.migration_speed[new_id] = parameters.motility_speed;
	data.migration_bias[new_id] = parameters.motility_bias;
	data.restrict_to_2d[new_id] = static_cast<std::uint8_t>(is_2D);

	const index_t substrate_count = data.substrates_count;
	index_t chosen_substrate = 0;
	bool found_substrate = false;

	for (index_t s = 0; s < substrate_count; ++s)
	{
		const bool enabled = s < parameters.chemotaxis_enabled.size() ? parameters.chemotaxis_enabled[s] : false;
		const real_t sensitivity =
			enabled && s < parameters.chemotaxis_sensitivity.size() ? parameters.chemotaxis_sensitivity[s] : 0.0;

		data.chemotactic_sensitivities[new_id * substrate_count + s] = sensitivity;

		if (!found_substrate && enabled)
		{
			chosen_substrate = s;
			found_substrate = true;
		}
	}

	data.chemotaxis_index[new_id] = found_substrate ? chosen_substrate : 0;
	data.chemotaxis_direction[new_id] = found_substrate ? 1 : 0;

	data.neighbors[new_id].clear();
	data.springs[new_id].clear();
	data.attached_cells[new_id].clear();

	data.simple_pressure[new_id] = 0.0;

	return agent;
}

} // namespace physicore::mechanics::physicell
