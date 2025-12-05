#include "physicell/mechanical_agent.h"

#include <algorithm>
#include <cassert>

namespace physicore::mechanics::physicell {

namespace {
template <typename T>
T value_or(const std::vector<T>& values, index_t idx, T fallback)
{
	return idx < values.size() ? static_cast<T>(values[idx]) : fallback;
}
} // namespace

std::unique_ptr<mechanical_agent> mechanical_agent::add(index_t new_id, index_t cell_type,
														SimulationParameters& parameters)
{
	return add(new_id, cell_type, parameters, data);
}

std::unique_ptr<mechanical_agent> mechanical_agent::add(index_t new_id, index_t cell_type,
														SimulationParameters& parameters, agent_data& data)
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

	data.cell_cell_adhesion_strength[new_id] =
		value_or(parameters.cell_cell_adhesion_strength, cell_type, 0.0);
	data.cell_BM_adhesion_strength[new_id] = value_or(parameters.cell_BM_adhesion_strength, cell_type, 0.0);
	data.cell_cell_repulsion_strength[new_id] =
		value_or(parameters.cell_cell_repulsion_strength, cell_type, 0.0);
	data.cell_BM_repulsion_strength[new_id] = value_or(parameters.cell_BM_repulsion_strength, cell_type, 0.0);

	const index_t type_count = data.agent_types_count;
	for (index_t other = 0; other < type_count; ++other)
	{
		const index_t param_idx = cell_type * type_count + other;
		const real_t affinity =
			param_idx < parameters.cell_adhesion_affinity.size() ? parameters.cell_adhesion_affinity[param_idx] : 0.0;
		data.cell_adhesion_affinities[new_id * type_count + other] = affinity;
	}

	data.relative_maximum_adhesion_distance[new_id] =
		value_or(parameters.relative_maximum_adhesion_distance, cell_type, 0.0);
	data.maximum_number_of_attachments[new_id] =
		value_or(parameters.maximum_number_of_attachments, cell_type, index_t(0));
	data.attachment_elastic_constant[new_id] =
		value_or(parameters.attachment_elastic_coefficient, cell_type, 0.0);
	data.attachment_rate[new_id] = value_or(parameters.attachment_rate, cell_type, 0.0);
	data.detachment_rate[new_id] = value_or(parameters.detachment_rate, cell_type, 0.0);

	const bool motile = value_or(parameters.is_movable, cell_type, false);
	data.is_motile[new_id] = static_cast<std::uint8_t>(motile);
	data.is_movable[new_id] = static_cast<std::uint8_t>(motile);
	data.persistence_time[new_id] = value_or(parameters.motility_persistence_time, cell_type, 0.0);
	data.migration_speed[new_id] = value_or(parameters.motility_speed, cell_type, 0.0);
	data.migration_bias[new_id] = value_or(parameters.motility_bias, cell_type, 0.0);
	data.restrict_to_2d[new_id] = static_cast<std::uint8_t>(parameters.is_2D);

	const index_t substrate_count = data.substrates_count;
	index_t chosen_substrate = 0;
	bool found_substrate = false;

	for (index_t s = 0; s < substrate_count; ++s)
	{
		const index_t param_idx = cell_type * substrate_count + s;
		const bool enabled =
			param_idx < parameters.chemotaxis_enabled.size() ? parameters.chemotaxis_enabled[param_idx] : false;
		const real_t sensitivity =
			enabled && param_idx < parameters.chemotaxis_sensitivity.size()
				? parameters.chemotaxis_sensitivity[param_idx]
				: 0.0;

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
