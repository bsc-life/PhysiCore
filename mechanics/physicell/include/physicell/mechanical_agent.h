#pragma once

#include <memory>
#include <span>
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
	agent_data& data;

public:
	using DataType = agent_data;
	using InterfaceType = agent_interface;

	mechanical_agent(index_t index, agent_data& data)
		: base_agent_interface(index),
		  physicore::base_agent_generic_storage<base_agent_data>(index, data.base_data),
		  data(data)
	{}

	mechanical_agent(index_t index, std::tuple<std::unique_ptr<base_agent_data>, std::unique_ptr<agent_data>>& datas)
		: mechanical_agent(index, *std::get<std::unique_ptr<agent_data>>(datas))
	{}

	std::span<real_t> velocity() override
	{
		const index_t dims = data.base_data.dims;
		return std::span<real_t>(&data.velocity[this->index * dims], dims);
	}

	std::span<real_t> previous_velocity() override
	{
		const index_t dims = data.base_data.dims;
		return std::span<real_t>(&data.previous_velocity[this->index * dims], dims);
	}

	real_t& radius() override { return data.radius[this->index]; }

	real_t& cell_cell_adhesion_strength() override { return data.mechanics.cell_cell_adhesion_strength[this->index]; }
	real_t& cell_BM_adhesion_strength() override { return data.mechanics.cell_BM_adhesion_strength[this->index]; }
	real_t& cell_cell_repulsion_strength() override { return data.mechanics.cell_cell_repulsion_strength[this->index]; }
	real_t& cell_BM_repulsion_strength() override { return data.mechanics.cell_BM_repulsion_strength[this->index]; }

	std::span<real_t> cell_adhesion_affinities() override
	{
		return std::span<real_t>(&data.mechanics.cell_adhesion_affinities[this->index * data.agent_types_count],
								 data.agent_types_count);
	}

	real_t& relative_maximum_adhesion_distance() override
	{
		return data.mechanics.relative_maximum_adhesion_distance[this->index];
	}

	index_t& maximum_number_of_attachments() override
	{
		return data.mechanics.maximum_number_of_attachments[this->index];
	}
	real_t& attachment_elastic_constant() override { return data.mechanics.attachment_elastic_constant[this->index]; }
	real_t& attachment_rate() override { return data.mechanics.attachment_rate[this->index]; }
	real_t& detachment_rate() override { return data.mechanics.detachment_rate[this->index]; }

	std::uint8_t& is_motile() override { return data.motility.is_motile[this->index]; }
	real_t& persistence_time() override { return data.motility.persistence_time[this->index]; }
	real_t& migration_speed() override { return data.motility.migration_speed[this->index]; }

	std::span<real_t> migration_bias_direction() override
	{
		const index_t dims = data.base_data.dims;
		return std::span<real_t>(&data.motility.migration_bias_direction[this->index * dims], dims);
	}

	real_t& migration_bias() override { return data.motility.migration_bias[this->index]; }

	std::span<real_t> motility_vector() override
	{
		const index_t dims = data.base_data.dims;
		return std::span<real_t>(&data.motility.motility_vector[this->index * dims], dims);
	}

	std::uint8_t& restrict_to_2d() override { return data.motility.restrict_to_2d[this->index]; }
	index_t& chemotaxis_index() override { return data.motility.chemotaxis_index[this->index]; }
	index_t& chemotaxis_direction() override { return data.motility.chemotaxis_direction[this->index]; }

	std::span<real_t> chemotactic_sensitivities() override
	{
		return std::span<real_t>(&data.motility.chemotactic_sensitivities[this->index * data.substrates_count],
								 data.substrates_count);
	}

	std::vector<index_t>& neighbors() override { return data.state.neighbors[this->index]; }
	std::vector<index_t>& springs() override { return data.state.springs[this->index]; }
	std::vector<index_t>& attached_cells() override { return data.state.attached_cells[this->index]; }

	std::span<real_t> orientation() override
	{
		const index_t dims = data.base_data.dims;
		return std::span<real_t>(&data.state.orientation[this->index * dims], dims);
	}

	real_t& simple_pressure() override { return data.state.simple_pressure[this->index]; }
	index_t& agent_type_index() override { return data.state.agent_type_index[this->index]; }
	std::uint8_t& is_movable() override { return data.state.is_movable[this->index]; }

};

} // namespace physicore::mechanics::physicell
