#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include <common/base_agent_interface.h>
#include <common/types.h>

namespace physicore::mechanics::physicell {

class agent_interface : public virtual base_agent_interface
{
public:
	virtual std::span<real_t> velocity() = 0;
	virtual std::span<real_t> previous_velocity() = 0;
	virtual real_t& radius() = 0;

	virtual real_t& cell_cell_adhesion_strength() = 0;
	virtual real_t& cell_BM_adhesion_strength() = 0;
	virtual real_t& cell_cell_repulsion_strength() = 0;
	virtual real_t& cell_BM_repulsion_strength() = 0;

	virtual std::span<real_t> cell_adhesion_affinities() = 0;
	virtual real_t& relative_maximum_adhesion_distance() = 0;
	virtual index_t& maximum_number_of_attachments() = 0;
	virtual real_t& attachment_elastic_constant() = 0;
	virtual real_t& attachment_rate() = 0;
	virtual real_t& detachment_rate() = 0;

	virtual std::uint8_t& is_motile() = 0;
	virtual real_t& persistence_time() = 0;
	virtual real_t& migration_speed() = 0;
	virtual std::span<real_t> migration_bias_direction() = 0;
	virtual real_t& migration_bias() = 0;
	virtual std::span<real_t> motility_vector() = 0;
	virtual std::uint8_t& restrict_to_2d() = 0;

	virtual index_t& chemotaxis_index() = 0;
	virtual index_t& chemotaxis_direction() = 0;
	virtual std::span<real_t> chemotactic_sensitivities() = 0;

	virtual std::vector<index_t>& neighbors() = 0;
	virtual std::vector<index_t>& springs() = 0;
	virtual std::vector<index_t>& attached_cells() = 0;

	virtual std::span<real_t> orientation() = 0;
	virtual real_t& simple_pressure() = 0;
	virtual index_t& agent_type_index() = 0;
	virtual std::uint8_t& is_movable() = 0;
};

} // namespace physicore::mechanics::physicell
