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
	// Sync with base class
	agents_count = base_data.agents_count;

	resize_storage();
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
		motility_data.direction_update_funcs[position] = std::move(motility_data.direction_update_funcs[last]);

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

	// Sync with base class
	agents_count = base_data.agents_count;

	resize_storage();
}

void mechanical_agent_data::resize_storage()
{
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
	motility_data.direction_update_funcs.resize(agents_count);
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

} // namespace physicore::mechanics::physicell
