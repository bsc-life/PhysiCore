#include "physicell/mechanical_agent.h"

#include <tuple>

namespace physicore::mechanics::physicell {

mechanical_agent::mechanical_agent(index_t index, mechanical_agent_data& data)
	: base_agent_interface(index),
	  physicore::base_agent_generic_storage<base_agent_data>(index, data.base_data),
	  data(data)
{}

mechanical_agent::mechanical_agent(
	index_t index, std::tuple<std::unique_ptr<base_agent_data>, std::unique_ptr<mechanical_agent_data>>& datas)
	: mechanical_agent(index, *std::get<std::unique_ptr<mechanical_agent_data>>(datas))
{}

std::span<real_t> mechanical_agent::velocity()
{
	const index_t dims = data.base_data.dims;
	return { &data.velocity[this->index * dims], dims };
}

std::span<real_t> mechanical_agent::previous_velocity()
{
	const index_t dims = data.base_data.dims;
	return { &data.previous_velocity[this->index * dims], dims };
}

real_t& mechanical_agent::radius() { return data.radius[this->index]; }

real_t& mechanical_agent::cell_cell_adhesion_strength()
{
	return data.mechanics_data.cell_cell_adhesion_strength[this->index];
}

real_t& mechanical_agent::cell_BM_adhesion_strength()
{
	return data.mechanics_data.cell_BM_adhesion_strength[this->index];
}

real_t& mechanical_agent::cell_cell_repulsion_strength()
{
	return data.mechanics_data.cell_cell_repulsion_strength[this->index];
}

real_t& mechanical_agent::cell_BM_repulsion_strength()
{
	return data.mechanics_data.cell_BM_repulsion_strength[this->index];
}

std::span<real_t> mechanical_agent::cell_adhesion_affinities()
{
	return { &data.mechanics_data.cell_adhesion_affinities[this->index * data.agent_types_count],
			 data.agent_types_count };
}

real_t& mechanical_agent::relative_maximum_adhesion_distance()
{
	return data.mechanics_data.relative_maximum_adhesion_distance[this->index];
}

index_t& mechanical_agent::maximum_number_of_attachments()
{
	return data.mechanics_data.maximum_number_of_attachments[this->index];
}

real_t& mechanical_agent::attachment_elastic_constant()
{
	return data.mechanics_data.attachment_elastic_constant[this->index];
}

real_t& mechanical_agent::attachment_rate() { return data.mechanics_data.attachment_rate[this->index]; }

real_t& mechanical_agent::detachment_rate() { return data.mechanics_data.detachment_rate[this->index]; }

std::uint8_t& mechanical_agent::is_motile() { return data.motility_data.is_motile[this->index]; }

real_t& mechanical_agent::persistence_time() { return data.motility_data.persistence_time[this->index]; }

real_t& mechanical_agent::migration_speed() { return data.motility_data.migration_speed[this->index]; }

std::span<real_t> mechanical_agent::migration_bias_direction()
{
	const index_t dims = data.base_data.dims;
	return { &data.motility_data.migration_bias_direction[this->index * dims], dims };
}

real_t& mechanical_agent::migration_bias() { return data.motility_data.migration_bias[this->index]; }

std::span<real_t> mechanical_agent::motility_vector()
{
	const index_t dims = data.base_data.dims;
	return { &data.motility_data.motility_vector[this->index * dims], dims };
}

std::uint8_t& mechanical_agent::restrict_to_2d() { return data.motility_data.restrict_to_2d[this->index]; }

index_t& mechanical_agent::chemotaxis_index() { return data.motility_data.chemotaxis_index[this->index]; }

index_t& mechanical_agent::chemotaxis_direction() { return data.motility_data.chemotaxis_direction[this->index]; }

std::span<real_t> mechanical_agent::chemotactic_sensitivities()
{
	return { &data.motility_data.chemotactic_sensitivities[this->index * data.substrates_count],
			 data.substrates_count };
}

std::span<index_t> mechanical_agent::neighbors() { return std::span<index_t>(data.state_data.neighbors[this->index]); }

std::span<index_t> mechanical_agent::springs() { return std::span<index_t>(data.state_data.springs[this->index]); }

std::span<index_t> mechanical_agent::attached_cells()
{
	return std::span<index_t>(data.state_data.attached_cells[this->index]);
}

std::span<real_t> mechanical_agent::orientation()
{
	const index_t dims = data.base_data.dims;
	return { &data.state_data.orientation[this->index * dims], dims };
}

real_t& mechanical_agent::simple_pressure() { return data.state_data.simple_pressure[this->index]; }

index_t& mechanical_agent::agent_type_index() { return data.state_data.agent_type_index[this->index]; }

std::uint8_t& mechanical_agent::is_movable() { return data.state_data.is_movable[this->index]; }

} // namespace physicore::mechanics::physicell
