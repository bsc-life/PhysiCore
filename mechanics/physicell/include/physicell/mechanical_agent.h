#pragma once

#include <memory>
#include <tuple>
#include <vector>

#include <common/base_agent_data.h>
#include <common/base_agent_generic_storage.h>
#include <common/types.h>

#include "agent_data.h"
#include "agent_interface.h"
#include "mechanical_parameters.h"

namespace physicore::mechanics::physicell {

class mechanical_agent : public physicore::base_agent_generic_storage<base_agent_data>, public virtual agent_interface
{
protected:
	index_t index;
	index_t type_index;
	agent_data& data;

public:
	using DataType = agent_data;
	using InterfaceType = agent_interface;

	mechanical_agent(index_t index, agent_data& data)
		: base_agent_interface(index),
		  physicore::base_agent_generic_storage<base_agent_data>(index, data.base_data),
		  index(index),
		  type_index(0),
		  data(data)
	{}

	mechanical_agent(index_t index, std::tuple<std::unique_ptr<base_agent_data>, std::unique_ptr<agent_data>>& datas)
		: mechanical_agent(index, *std::get<std::unique_ptr<agent_data>>(datas))
	{}

	std::span<real_t> velocity() override
	{
		return std::span<real_t>(&data.velocities[this->index * data.dims], data.dims);
	}

	std::span<real_t> previous_velocity() override
	{
		return std::span<real_t>(&data.previous_velocities[this->index * data.dims], data.dims);
	}

	real_t& radius() override { return data.radius[this->index]; }

	real_t& cell_cell_adhesion_strength() override { return data.cell_cell_adhesion_strength[this->index]; }
	real_t& cell_BM_adhesion_strength() override { return data.cell_BM_adhesion_strength[this->index]; }
	real_t& cell_cell_repulsion_strength() override { return data.cell_cell_repulsion_strength[this->index]; }
	real_t& cell_BM_repulsion_strength() override { return data.cell_BM_repulsion_strength[this->index]; }

	std::span<real_t> cell_adhesion_affinities() override
	{
		return std::span<real_t>(&data.cell_adhesion_affinities[this->index * data.agent_types_count],
								 data.agent_types_count);
	}

	real_t& relative_maximum_adhesion_distance() override
	{
		return data.relative_maximum_adhesion_distance[this->index];
	}

	index_t& maximum_number_of_attachments() override { return data.maximum_number_of_attachments[this->index]; }
	real_t& attachment_elastic_constant() override { return data.attachment_elastic_constant[this->index]; }
	real_t& attachment_rate() override { return data.attachment_rate[this->index]; }
	real_t& detachment_rate() override { return data.detachment_rate[this->index]; }

	std::uint8_t& is_motile() override { return data.is_motile[this->index]; }
	real_t& persistence_time() override { return data.persistence_time[this->index]; }
	real_t& migration_speed() override { return data.migration_speed[this->index]; }

	std::span<real_t> migration_bias_direction() override
	{
		return std::span<real_t>(&data.migration_bias_direction[this->index * data.dims], data.dims);
	}

	real_t& migration_bias() override { return data.migration_bias[this->index]; }

	std::span<real_t> motility_vector() override
	{
		return std::span<real_t>(&data.motility_vector[this->index * data.dims], data.dims);
	}

	std::uint8_t& restrict_to_2d() override { return data.restrict_to_2d[this->index]; }
	index_t& chemotaxis_index() override { return data.chemotaxis_index[this->index]; }
	index_t& chemotaxis_direction() override { return data.chemotaxis_direction[this->index]; }

	std::span<real_t> chemotactic_sensitivities() override
	{
		return std::span<real_t>(&data.chemotactic_sensitivities[this->index * data.substrates_count],
								 data.substrates_count);
	}

	std::vector<index_t>& neighbors() override { return data.neighbors[this->index]; }
	std::vector<index_t>& springs() override { return data.springs[this->index]; }
	std::vector<index_t>& attached_cells() override { return data.attached_cells[this->index]; }

	std::span<real_t> orientation() override
	{
		return std::span<real_t>(&data.orientation[this->index * data.dims], data.dims);
	}

	real_t& simple_pressure() override { return data.simple_pressure[this->index]; }
	index_t& cell_definition_index() override { return data.cell_definition_indices[this->index]; }
	std::uint8_t& is_movable() override { return data.is_movable[this->index]; }

	std::unique_ptr<mechanical_agent> add(index_t new_id, index_t cell_type, mechanical_parameters& parameters,
										  bool is_2D);
	static std::unique_ptr<mechanical_agent> add(index_t new_id, index_t cell_type, mechanical_parameters& parameters,
												 agent_data& data, bool is_2D);
};

} // namespace physicore::mechanics::physicell
