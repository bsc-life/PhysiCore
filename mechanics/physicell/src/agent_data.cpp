#include "physicell/agent_data.h"

#include <algorithm>
#include <cassert>
#include <memory>

#include "physicell/mechanical_agent.h"
#include "physicell/mechanical_parameters.h"

namespace physicore::mechanics::physicell {

mechanical_agent_data::mechanical_agent_data(physicore::base_agent_data_generic_storage<std::vector>& base_data,
											 index_t agent_types_count, index_t substrates_count)
	: base_data(base_data), agent_types_count(agent_types_count), substrates_count(substrates_count)
{}

void mechanical_agent_data::add()
{
	agents_count = base_data.agents_count;

	velocity.resize(agents_count * base_data.dims);
	previous_velocity.resize(agents_count * base_data.dims);
	radius.resize(agents_count);

	// Mechanics properties
	mechanics.cell_cell_adhesion_strength.resize(agents_count);
	mechanics.cell_BM_adhesion_strength.resize(agents_count);
	mechanics.cell_cell_repulsion_strength.resize(agents_count);
	mechanics.cell_BM_repulsion_strength.resize(agents_count);
	mechanics.cell_adhesion_affinities.resize(agents_count * agent_types_count);
	mechanics.relative_maximum_adhesion_distance.resize(agents_count);
	mechanics.maximum_number_of_attachments.resize(agents_count);
	mechanics.attachment_elastic_constant.resize(agents_count);
	mechanics.attachment_rate.resize(agents_count);
	mechanics.detachment_rate.resize(agents_count);

	// Motility properties
	motility.is_motile.resize(agents_count);
	motility.persistence_time.resize(agents_count);
	motility.migration_speed.resize(agents_count);
	motility.migration_bias_direction.resize(agents_count * base_data.dims);
	motility.migration_bias.resize(agents_count);
	motility.motility_vector.resize(agents_count * base_data.dims);
	motility.restrict_to_2d.resize(agents_count);
	motility.chemotaxis_index.resize(agents_count);
	motility.chemotaxis_direction.resize(agents_count);
	motility.chemotactic_sensitivities.resize(agents_count * substrates_count);

	// State properties
	state.neighbors.resize(agents_count);
	state.springs.resize(agents_count);
	state.attached_cells.resize(agents_count);
	state.orientation.resize(agents_count * base_data.dims);
	state.simple_pressure.resize(agents_count);
	state.agent_type_index.resize(agents_count);
	state.is_movable.resize(agents_count);
}

void mechanical_agent_data::remove_at(index_t position)
{
	assert(position < agents_count);
	if (position >= agents_count)
		return;

	using base_storage_t = physicore::base_agent_data_generic_storage<std::vector>;

	const index_t last = agents_count - 1;

	if (position != last)
	{
		base_storage_t::move_vector(&velocity[position * base_data.dims], &velocity[last * base_data.dims],
									base_data.dims);
		base_storage_t::move_vector(&previous_velocity[position * base_data.dims],
									&previous_velocity[last * base_data.dims], base_data.dims);
		base_storage_t::move_scalar(&radius[position], &radius[last]);

		// Mechanics properties
		base_storage_t::move_scalar(&mechanics.cell_cell_adhesion_strength[position],
									&mechanics.cell_cell_adhesion_strength[last]);
		base_storage_t::move_scalar(&mechanics.cell_BM_adhesion_strength[position],
									&mechanics.cell_BM_adhesion_strength[last]);
		base_storage_t::move_scalar(&mechanics.cell_cell_repulsion_strength[position],
									&mechanics.cell_cell_repulsion_strength[last]);
		base_storage_t::move_scalar(&mechanics.cell_BM_repulsion_strength[position],
									&mechanics.cell_BM_repulsion_strength[last]);

		base_storage_t::move_vector(&mechanics.cell_adhesion_affinities[position * agent_types_count],
									&mechanics.cell_adhesion_affinities[last * agent_types_count], agent_types_count);

		base_storage_t::move_scalar(&mechanics.relative_maximum_adhesion_distance[position],
									&mechanics.relative_maximum_adhesion_distance[last]);
		base_storage_t::move_scalar(&mechanics.maximum_number_of_attachments[position],
									&mechanics.maximum_number_of_attachments[last]);
		base_storage_t::move_scalar(&mechanics.attachment_elastic_constant[position],
									&mechanics.attachment_elastic_constant[last]);
		base_storage_t::move_scalar(&mechanics.attachment_rate[position], &mechanics.attachment_rate[last]);
		base_storage_t::move_scalar(&mechanics.detachment_rate[position], &mechanics.detachment_rate[last]);

		// Motility properties
		base_storage_t::move_scalar(&motility.is_motile[position], &motility.is_motile[last]);
		base_storage_t::move_scalar(&motility.persistence_time[position], &motility.persistence_time[last]);
		base_storage_t::move_scalar(&motility.migration_speed[position], &motility.migration_speed[last]);
		base_storage_t::move_vector(&motility.migration_bias_direction[position * base_data.dims],
									&motility.migration_bias_direction[last * base_data.dims], base_data.dims);
		base_storage_t::move_scalar(&motility.migration_bias[position], &motility.migration_bias[last]);
		base_storage_t::move_vector(&motility.motility_vector[position * base_data.dims],
									&motility.motility_vector[last * base_data.dims], base_data.dims);
		base_storage_t::move_scalar(&motility.restrict_to_2d[position], &motility.restrict_to_2d[last]);
		base_storage_t::move_scalar(&motility.chemotaxis_index[position], &motility.chemotaxis_index[last]);
		base_storage_t::move_scalar(&motility.chemotaxis_direction[position], &motility.chemotaxis_direction[last]);
		base_storage_t::move_vector(&motility.chemotactic_sensitivities[position * substrates_count],
									&motility.chemotactic_sensitivities[last * substrates_count], substrates_count);

		// State properties
		state.neighbors[position] = std::move(state.neighbors[last]);
		state.springs[position] = std::move(state.springs[last]);
		state.attached_cells[position] = std::move(state.attached_cells[last]);
		base_storage_t::move_vector(&state.orientation[position * base_data.dims],
									&state.orientation[last * base_data.dims], base_data.dims);
		base_storage_t::move_scalar(&state.simple_pressure[position], &state.simple_pressure[last]);
		base_storage_t::move_scalar(&state.agent_type_index[position], &state.agent_type_index[last]);
		base_storage_t::move_scalar(&state.is_movable[position], &state.is_movable[last]);
	}

	agents_count = base_data.agents_count;

	velocity.resize(agents_count * base_data.dims);
	previous_velocity.resize(agents_count * base_data.dims);
	radius.resize(agents_count);

	// Mechanics properties
	mechanics.cell_cell_adhesion_strength.resize(agents_count);
	mechanics.cell_BM_adhesion_strength.resize(agents_count);
	mechanics.cell_cell_repulsion_strength.resize(agents_count);
	mechanics.cell_BM_repulsion_strength.resize(agents_count);
	mechanics.cell_adhesion_affinities.resize(agents_count * agent_types_count);
	mechanics.relative_maximum_adhesion_distance.resize(agents_count);
	mechanics.maximum_number_of_attachments.resize(agents_count);
	mechanics.attachment_elastic_constant.resize(agents_count);
	mechanics.attachment_rate.resize(agents_count);
	mechanics.detachment_rate.resize(agents_count);

	// Motility properties
	motility.is_motile.resize(agents_count);
	motility.persistence_time.resize(agents_count);
	motility.migration_speed.resize(agents_count);
	motility.migration_bias_direction.resize(agents_count * base_data.dims);
	motility.migration_bias.resize(agents_count);
	motility.motility_vector.resize(agents_count * base_data.dims);
	motility.restrict_to_2d.resize(agents_count);
	motility.chemotaxis_index.resize(agents_count);
	motility.chemotaxis_direction.resize(agents_count);
	motility.chemotactic_sensitivities.resize(agents_count * substrates_count);

	// State properties
	state.neighbors.resize(agents_count);
	state.springs.resize(agents_count);
	state.attached_cells.resize(agents_count);
	state.orientation.resize(agents_count * base_data.dims);
	state.simple_pressure.resize(agents_count);
	state.agent_type_index.resize(agents_count);
	state.is_movable.resize(agents_count);
}



void mechanical_agent_data::add(index_t new_id, index_t cell_type, mechanical_parameters& parameters, bool is_2D)
{
	assert(new_id == base_data.agents_count);

	base_data.add();
	add();

	assert(new_id < agents_count);

	const index_t dim_offset = new_id * base_data.dims;
	std::fill_n(&velocity[dim_offset], base_data.dims, 0.0);
	std::fill_n(&previous_velocity[dim_offset], base_data.dims, 0.0);
	std::fill_n(&motility.migration_bias_direction[dim_offset], base_data.dims, 0.0);
	std::fill_n(&motility.motility_vector[dim_offset], base_data.dims, 0.0);
	std::fill_n(&state.orientation[dim_offset], base_data.dims, 0.0);

	radius[new_id] = 0.0;

	state.agent_type_index[new_id] = cell_type;

	// Directly assign scalar values from mechanical_parameters
	mechanics.cell_cell_adhesion_strength[new_id] = parameters.cell_cell_adhesion_strength;
	mechanics.cell_BM_adhesion_strength[new_id] = parameters.cell_BM_adhesion_strength;
	mechanics.cell_cell_repulsion_strength[new_id] = parameters.cell_cell_repulsion_strength;
	mechanics.cell_BM_repulsion_strength[new_id] = parameters.cell_BM_repulsion_strength;

	// Copy cell_adhesion_affinities vector directly
	const index_t affinity_offset = new_id * agent_types_count;
	index_t other = 0;
	for (; other < agent_types_count && other < static_cast<index_t>(parameters.cell_adhesion_affinity.size()); ++other)
	{
		mechanics.cell_adhesion_affinities[affinity_offset + other] = parameters.cell_adhesion_affinity[other];
	}
	for (; other < agent_types_count; ++other)
	{
		mechanics.cell_adhesion_affinities[affinity_offset + other] = 0.0;
	}

	mechanics.relative_maximum_adhesion_distance[new_id] = parameters.relative_maximum_adhesion_distance;
	mechanics.maximum_number_of_attachments[new_id] = parameters.maximum_number_of_attachments;
	mechanics.attachment_elastic_constant[new_id] = parameters.attachment_elastic_coefficient;
	mechanics.attachment_rate[new_id] = parameters.attachment_rate;
	mechanics.detachment_rate[new_id] = parameters.detachment_rate;

	motility.is_motile[new_id] = static_cast<std::uint8_t>(parameters.is_movable);
	state.is_movable[new_id] = static_cast<std::uint8_t>(parameters.is_movable);
	motility.persistence_time[new_id] = parameters.motility_persistence_time;
	motility.migration_speed[new_id] = parameters.motility_speed;
	motility.migration_bias[new_id] = parameters.motility_bias;
	motility.restrict_to_2d[new_id] = static_cast<std::uint8_t>(is_2D);

	index_t chosen_substrate = 0;
	bool found_substrate = false;

	for (index_t s = 0; s < substrates_count; ++s)
	{
		const bool enabled = s < parameters.chemotaxis_enabled.size() ? parameters.chemotaxis_enabled[s] : false;
		const real_t sensitivity =
			enabled && s < parameters.chemotaxis_sensitivity.size() ? parameters.chemotaxis_sensitivity[s] : 0.0;

		motility.chemotactic_sensitivities[new_id * substrates_count + s] = sensitivity;

		if (!found_substrate && enabled)
		{
			chosen_substrate = s;
			found_substrate = true;
		}
	}

	motility.chemotaxis_index[new_id] = found_substrate ? chosen_substrate : 0;
	motility.chemotaxis_direction[new_id] = found_substrate ? 1 : 0;

	state.neighbors[new_id].clear();
	state.springs[new_id].clear();
	state.attached_cells[new_id].clear();

	state.simple_pressure[new_id] = 0.0;
}

} // namespace physicore::mechanics::physicell
