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

class mechanical_agent : public physicore::base_agent_generic_storage<base_agent_data>,
						 public virtual mechanical_agent_interface
{
protected:
	mechanical_agent_data& data;

public:
	using DataType = mechanical_agent_data;
	using InterfaceType = mechanical_agent_interface;

	mechanical_agent(index_t index, mechanical_agent_data& data);
	mechanical_agent(index_t index,
					 std::tuple<std::unique_ptr<base_agent_data>, std::unique_ptr<mechanical_agent_data>>& datas);

	std::span<real_t> velocity() override;
	std::span<real_t> previous_velocity() override;

	real_t& radius() override;

	real_t& cell_cell_adhesion_strength() override;
	real_t& cell_BM_adhesion_strength() override;
	real_t& cell_cell_repulsion_strength() override;
	real_t& cell_BM_repulsion_strength() override;

	std::span<real_t> cell_adhesion_affinities() override;
	real_t& relative_maximum_adhesion_distance() override;

	index_t& maximum_number_of_attachments() override;
	real_t& attachment_elastic_constant() override;
	real_t& attachment_rate() override;
	real_t& detachment_rate() override;

	std::uint8_t& is_motile() override;
	real_t& persistence_time() override;
	real_t& migration_speed() override;

	std::span<real_t> migration_bias_direction() override;
	real_t& migration_bias() override;
	std::span<real_t> motility_vector() override;

	std::uint8_t& restrict_to_2d() override;
	index_t& chemotaxis_index() override;
	index_t& chemotaxis_direction() override;

	std::span<real_t> chemotactic_sensitivities() override;

	std::span<index_t> neighbors() override;
	std::span<index_t> springs() override;
	std::span<index_t> attached_cells() override;

	std::span<real_t> orientation() override;

	real_t& simple_pressure() override;
	index_t& agent_type_index() override;
	std::uint8_t& is_movable() override;
};

} // namespace physicore::mechanics::physicell
