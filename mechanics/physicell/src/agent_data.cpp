#include "physicell/agent_data.h"

namespace physicore::mechanics::physicell {


template <template <typename...> typename ContainerType>
mechanical_agent_data<ContainerType>::mechanical_agent_data(
	physicore::base_agent_data_generic_storage<ContainerType>& base_data, index_t dims, index_t agent_types_count,
	index_t substrates_count)
	: base_data(base_data),
	  dims(base_data.dims <= 0 ? dims : base_data.dims),
	  agent_types_count(agent_types_count),
	  substrates_count(substrates_count)
{}

template <template <typename...> typename ContainerType>
void mechanical_agent_data<ContainerType>::add()
{
	agents_count = base_data.agents_count;
	dims = base_data.dims;

	velocities.resize(agents_count * dims);
	previous_velocities.resize(agents_count * dims);
	radius.resize(agents_count);

	cell_cell_adhesion_strength.resize(agents_count);
	cell_BM_adhesion_strength.resize(agents_count);
	cell_cell_repulsion_strength.resize(agents_count);
	cell_BM_repulsion_strength.resize(agents_count);
	cell_adhesion_affinities.resize(agents_count * agent_types_count);
	relative_maximum_adhesion_distance.resize(agents_count);
	maximum_number_of_attachments.resize(agents_count);
	attachment_elastic_constant.resize(agents_count);
	attachment_rate.resize(agents_count);
	detachment_rate.resize(agents_count);


	is_motile.resize(agents_count);
	persistence_time.resize(agents_count);
	migration_speed.resize(agents_count);
	migration_bias_direction.resize(agents_count * dims);
	migration_bias.resize(agents_count);
	motility_vector.resize(agents_count * dims);
	restrict_to_2d.resize(agents_count);
	chemotaxis_index.resize(agents_count);
	chemotaxis_direction.resize(agents_count);
	chemotactic_sensitivities.resize(agents_count * substrates_count);

	neighbors.resize(agents_count);
	springs.resize(agents_count);
	attached_cells.resize(agents_count);
	orientation.resize(agents_count * dims);
	simple_pressure.resize(agents_count);
	cell_definition_indices.resize(agents_count);
	is_movable.resize(agents_count);
}

template <template <typename...> typename ContainerType>
void mechanical_agent_data<ContainerType>::remove_at(index_t position)
{
	assert(position < agents_count);
	if (position >= agents_count)
		return;

	using base_storage_t = physicore::base_agent_data_generic_storage<ContainerType>;

	const index_t last = agents_count - 1;

	if (position != last)
	{
		base_storage_t::move_vector(&velocities[position * dims], &velocities[last * dims], dims);
		base_storage_t::move_vector(&previous_velocities[position * dims], &previous_velocities[last * dims], dims);
		base_storage_t::move_scalar(&radius[position], &radius[last]);

		base_storage_t::move_scalar(&cell_cell_adhesion_strength[position], &cell_cell_adhesion_strength[last]);
		base_storage_t::move_scalar(&cell_BM_adhesion_strength[position], &cell_BM_adhesion_strength[last]);
		base_storage_t::move_scalar(&cell_cell_repulsion_strength[position], &cell_cell_repulsion_strength[last]);
		base_storage_t::move_scalar(&cell_BM_repulsion_strength[position], &cell_BM_repulsion_strength[last]);

		base_storage_t::move_vector(&cell_adhesion_affinities[position * agent_types_count],
									&cell_adhesion_affinities[last * agent_types_count], agent_types_count);

		base_storage_t::move_scalar(&relative_maximum_adhesion_distance[position],
									&relative_maximum_adhesion_distance[last]);
		base_storage_t::move_scalar(&maximum_number_of_attachments[position], &maximum_number_of_attachments[last]);
		base_storage_t::move_scalar(&attachment_elastic_constant[position], &attachment_elastic_constant[last]);
		base_storage_t::move_scalar(&attachment_rate[position], &attachment_rate[last]);
		base_storage_t::move_scalar(&detachment_rate[position], &detachment_rate[last]);

		base_storage_t::move_scalar(&is_motile[position], &is_motile[last]);
		base_storage_t::move_scalar(&persistence_time[position], &persistence_time[last]);
		base_storage_t::move_scalar(&migration_speed[position], &migration_speed[last]);
		base_storage_t::move_vector(&migration_bias_direction[position * dims], &migration_bias_direction[last * dims],
									dims);
		base_storage_t::move_scalar(&migration_bias[position], &migration_bias[last]);
		base_storage_t::move_vector(&motility_vector[position * dims], &motility_vector[last * dims], dims);
		base_storage_t::move_scalar(&restrict_to_2d[position], &restrict_to_2d[last]);
		base_storage_t::move_scalar(&chemotaxis_index[position], &chemotaxis_index[last]);
		base_storage_t::move_scalar(&chemotaxis_direction[position], &chemotaxis_direction[last]);
		base_storage_t::move_vector(&chemotactic_sensitivities[position * substrates_count],
									&chemotactic_sensitivities[last * substrates_count], substrates_count);

		neighbors[position] = std::move(neighbors[last]);
		springs[position] = std::move(springs[last]);
		attached_cells[position] = std::move(attached_cells[last]);
		base_storage_t::move_vector(&orientation[position * dims], &orientation[last * dims], dims);
		base_storage_t::move_scalar(&simple_pressure[position], &simple_pressure[last]);
		base_storage_t::move_scalar(&cell_definition_indices[position], &cell_definition_indices[last]);
		base_storage_t::move_scalar(&is_movable[position], &is_movable[last]);
	}

	agents_count = base_data.agents_count;

	velocities.resize(agents_count * dims);
	previous_velocities.resize(agents_count * dims);
	radius.resize(agents_count);

	cell_cell_adhesion_strength.resize(agents_count);
	cell_BM_adhesion_strength.resize(agents_count);
	cell_cell_repulsion_strength.resize(agents_count);
	cell_BM_repulsion_strength.resize(agents_count);
	cell_adhesion_affinities.resize(agents_count * agent_types_count);
	relative_maximum_adhesion_distance.resize(agents_count);
	maximum_number_of_attachments.resize(agents_count);
	attachment_elastic_constant.resize(agents_count);
	attachment_rate.resize(agents_count);
	detachment_rate.resize(agents_count);

	is_motile.resize(agents_count);
	persistence_time.resize(agents_count);
	migration_speed.resize(agents_count);
	migration_bias_direction.resize(agents_count * dims);
	migration_bias.resize(agents_count);
	motility_vector.resize(agents_count * dims);
	restrict_to_2d.resize(agents_count);
	chemotaxis_index.resize(agents_count);
	chemotaxis_direction.resize(agents_count);
	chemotactic_sensitivities.resize(agents_count * substrates_count);

	neighbors.resize(agents_count);
	springs.resize(agents_count);
	attached_cells.resize(agents_count);
	orientation.resize(agents_count * dims);
	simple_pressure.resize(agents_count);
	cell_definition_indices.resize(agents_count);
	is_movable.resize(agents_count);
}

template struct mechanical_agent_data<std::vector>;

} // namespace physicore::mechanics::physicell
