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
	mechanics_data.cell_cell_adhesion_strength.resize(agents_count);
	mechanics_data.cell_BM_adhesion_strength.resize(agents_count);
	mechanics_data.cell_cell_repulsion_strength.resize(agents_count);
	mechanics_data.cell_BM_repulsion_strength.resize(agents_count);
	mechanics_data.cell_adhesion_affinities.resize(agents_count * agent_types_count);
	mechanics_data.relative_maximum_adhesion_distance.resize(agents_count);
	mechanics_data.maximum_number_of_attachments.resize(agents_count);
	mechanics_data.attachment_elastic_constant.resize(agents_count);
	mechanics_data.attachment_rate.resize(agents_count);
	mechanics_data.detachment_rate.resize(agents_count);

	// Motility properties
	motility_data.is_motile.resize(agents_count);
	motility_data.persistence_time.resize(agents_count);
	motility_data.migration_speed.resize(agents_count);
	motility_data.migration_bias_direction.resize(agents_count * base_data.dims);
	motility_data.migration_bias.resize(agents_count);
	motility_data.motility_vector.resize(agents_count * base_data.dims);
	motility_data.restrict_to_2d.resize(agents_count);
	motility_data.chemotaxis_index.resize(agents_count);
	motility_data.chemotaxis_direction.resize(agents_count);
	motility_data.chemotactic_sensitivities.resize(agents_count * substrates_count);

	// State properties
	state_data.neighbors.resize(agents_count);
	state_data.springs.resize(agents_count);
	state_data.attached_cells.resize(agents_count);
	state_data.orientation.resize(agents_count * base_data.dims);
	state_data.simple_pressure.resize(agents_count);
	state_data.agent_type_index.resize(agents_count);
	state_data.is_movable.resize(agents_count);
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
		base_storage_t::move_scalar(&mechanics_data.cell_cell_adhesion_strength[position],
									&mechanics_data.cell_cell_adhesion_strength[last]);
		base_storage_t::move_scalar(&mechanics_data.cell_BM_adhesion_strength[position],
									&mechanics_data.cell_BM_adhesion_strength[last]);
		base_storage_t::move_scalar(&mechanics_data.cell_cell_repulsion_strength[position],
									&mechanics_data.cell_cell_repulsion_strength[last]);
		base_storage_t::move_scalar(&mechanics_data.cell_BM_repulsion_strength[position],
									&mechanics_data.cell_BM_repulsion_strength[last]);

		base_storage_t::move_vector(&mechanics_data.cell_adhesion_affinities[position * agent_types_count],
									&mechanics_data.cell_adhesion_affinities[last * agent_types_count],
									agent_types_count);

		base_storage_t::move_scalar(&mechanics_data.relative_maximum_adhesion_distance[position],
									&mechanics_data.relative_maximum_adhesion_distance[last]);
		base_storage_t::move_scalar(&mechanics_data.maximum_number_of_attachments[position],
									&mechanics_data.maximum_number_of_attachments[last]);
		base_storage_t::move_scalar(&mechanics_data.attachment_elastic_constant[position],
									&mechanics_data.attachment_elastic_constant[last]);
		base_storage_t::move_scalar(&mechanics_data.attachment_rate[position], &mechanics_data.attachment_rate[last]);
		base_storage_t::move_scalar(&mechanics_data.detachment_rate[position], &mechanics_data.detachment_rate[last]);

		// Motility properties
		base_storage_t::move_scalar(&motility_data.is_motile[position], &motility_data.is_motile[last]);
		base_storage_t::move_scalar(&motility_data.persistence_time[position], &motility_data.persistence_time[last]);
		base_storage_t::move_scalar(&motility_data.migration_speed[position], &motility_data.migration_speed[last]);
		base_storage_t::move_vector(&motility_data.migration_bias_direction[position * base_data.dims],
									&motility_data.migration_bias_direction[last * base_data.dims], base_data.dims);
		base_storage_t::move_scalar(&motility_data.migration_bias[position], &motility_data.migration_bias[last]);
		base_storage_t::move_vector(&motility_data.motility_vector[position * base_data.dims],
									&motility_data.motility_vector[last * base_data.dims], base_data.dims);
		base_storage_t::move_scalar(&motility_data.restrict_to_2d[position], &motility_data.restrict_to_2d[last]);
		base_storage_t::move_scalar(&motility_data.chemotaxis_index[position], &motility_data.chemotaxis_index[last]);
		base_storage_t::move_scalar(&motility_data.chemotaxis_direction[position],
									&motility_data.chemotaxis_direction[last]);
		base_storage_t::move_vector(&motility_data.chemotactic_sensitivities[position * substrates_count],
									&motility_data.chemotactic_sensitivities[last * substrates_count],
									substrates_count);

		// State properties
		state_data.neighbors[position] = std::move(state_data.neighbors[last]);
		state_data.springs[position] = std::move(state_data.springs[last]);
		state_data.attached_cells[position] = std::move(state_data.attached_cells[last]);
		base_storage_t::move_vector(&state_data.orientation[position * base_data.dims],
									&state_data.orientation[last * base_data.dims], base_data.dims);
		base_storage_t::move_scalar(&state_data.simple_pressure[position], &state_data.simple_pressure[last]);
		base_storage_t::move_scalar(&state_data.agent_type_index[position], &state_data.agent_type_index[last]);
		base_storage_t::move_scalar(&state_data.is_movable[position], &state_data.is_movable[last]);
	}

	agents_count = base_data.agents_count;

	velocity.resize(agents_count * base_data.dims);
	previous_velocity.resize(agents_count * base_data.dims);
	radius.resize(agents_count);

	// Mechanics properties
	mechanics_data.cell_cell_adhesion_strength.resize(agents_count);
	mechanics_data.cell_BM_adhesion_strength.resize(agents_count);
	mechanics_data.cell_cell_repulsion_strength.resize(agents_count);
	mechanics_data.cell_BM_repulsion_strength.resize(agents_count);
	mechanics_data.cell_adhesion_affinities.resize(agents_count * agent_types_count);
	mechanics_data.relative_maximum_adhesion_distance.resize(agents_count);
	mechanics_data.maximum_number_of_attachments.resize(agents_count);
	mechanics_data.attachment_elastic_constant.resize(agents_count);
	mechanics_data.attachment_rate.resize(agents_count);
	mechanics_data.detachment_rate.resize(agents_count);

	// Motility properties
	motility_data.is_motile.resize(agents_count);
	motility_data.persistence_time.resize(agents_count);
	motility_data.migration_speed.resize(agents_count);
	motility_data.migration_bias_direction.resize(agents_count * base_data.dims);
	motility_data.migration_bias.resize(agents_count);
	motility_data.motility_vector.resize(agents_count * base_data.dims);
	motility_data.restrict_to_2d.resize(agents_count);
	motility_data.chemotaxis_index.resize(agents_count);
	motility_data.chemotaxis_direction.resize(agents_count);
	motility_data.chemotactic_sensitivities.resize(agents_count * substrates_count);

	// State properties
	state_data.neighbors.resize(agents_count);
	state_data.springs.resize(agents_count);
	state_data.attached_cells.resize(agents_count);
	state_data.orientation.resize(agents_count * base_data.dims);
	state_data.simple_pressure.resize(agents_count);
	state_data.agent_type_index.resize(agents_count);
	state_data.is_movable.resize(agents_count);
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
	std::fill_n(&motility_data.migration_bias_direction[dim_offset], base_data.dims, 0.0);
	std::fill_n(&motility_data.motility_vector[dim_offset], base_data.dims, 0.0);
	std::fill_n(&state_data.orientation[dim_offset], base_data.dims, 0.0);

	radius[new_id] = 0.0;

	state_data.agent_type_index[new_id] = cell_type;

	// Directly assign scalar values from mechanical_parameters
	mechanics_data.cell_cell_adhesion_strength[new_id] = parameters.cell_cell_adhesion_strength;
	mechanics_data.cell_BM_adhesion_strength[new_id] = parameters.cell_BM_adhesion_strength;
	mechanics_data.cell_cell_repulsion_strength[new_id] = parameters.cell_cell_repulsion_strength;
	mechanics_data.cell_BM_repulsion_strength[new_id] = parameters.cell_BM_repulsion_strength;

	// Copy cell_adhesion_affinities vector directly
	const index_t affinity_offset = new_id * agent_types_count;
	index_t other = 0;
	for (; other < agent_types_count && other < static_cast<index_t>(parameters.cell_adhesion_affinity.size()); ++other)
	{
		mechanics_data.cell_adhesion_affinities[affinity_offset + other] = parameters.cell_adhesion_affinity[other];
	}
	for (; other < agent_types_count; ++other)
	{
		mechanics_data.cell_adhesion_affinities[affinity_offset + other] = 0.0;
	}

	mechanics_data.relative_maximum_adhesion_distance[new_id] = parameters.relative_maximum_adhesion_distance;
	mechanics_data.maximum_number_of_attachments[new_id] = parameters.maximum_number_of_attachments;
	mechanics_data.attachment_elastic_constant[new_id] = parameters.attachment_elastic_coefficient;
	mechanics_data.attachment_rate[new_id] = parameters.attachment_rate;
	mechanics_data.detachment_rate[new_id] = parameters.detachment_rate;

	motility_data.is_motile[new_id] = static_cast<std::uint8_t>(parameters.is_movable);
	state_data.is_movable[new_id] = static_cast<std::uint8_t>(parameters.is_movable);
	motility_data.persistence_time[new_id] = parameters.motility_persistence_time;
	motility_data.migration_speed[new_id] = parameters.motility_speed;
	motility_data.migration_bias[new_id] = parameters.motility_bias;
	motility_data.restrict_to_2d[new_id] = static_cast<std::uint8_t>(is_2D);

	index_t chosen_substrate = 0;
	bool found_substrate = false;

	for (index_t s = 0; s < substrates_count; ++s)
	{
		const bool enabled = s < parameters.chemotaxis_enabled.size() ? parameters.chemotaxis_enabled[s] : false;
		const real_t sensitivity =
			enabled && s < parameters.chemotaxis_sensitivity.size() ? parameters.chemotaxis_sensitivity[s] : 0.0;

		motility_data.chemotactic_sensitivities[new_id * substrates_count + s] = sensitivity;

		if (!found_substrate && enabled)
		{
			chosen_substrate = s;
			found_substrate = true;
		}
	}

	motility_data.chemotaxis_index[new_id] = found_substrate ? chosen_substrate : 0;
	motility_data.chemotaxis_direction[new_id] = found_substrate ? 1 : 0;

	state_data.neighbors[new_id].clear();
	state_data.springs[new_id].clear();
	state_data.attached_cells[new_id].clear();

	state_data.simple_pressure[new_id] = 0.0;
}

} // namespace physicore::mechanics::physicell
