#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <common/base_agent_data.h>
#include <common/types.h>

namespace physicore::mechanics::physicell {

template <template <typename...> typename ContainerType = std::vector>
struct mechanical_agent_data
{
public:
	physicore::base_agent_data_generic_storage<ContainerType>& base_data;

	index_t agents_count = 0;
	index_t dims = 3;
	index_t agent_types_count = 0;
	index_t substrates_count = 0;

	// Geometry / kinematics
	ContainerType<real_t> velocities;		   // current velocity per axis
	ContainerType<real_t> previous_velocities; // previous velocity per axis
	ContainerType<real_t> radius;			   // cell radius


	// Mechanics parameters (per agent; can be parsed from simulation parameters)
	ContainerType<real_t> cell_cell_adhesion_strength;
	ContainerType<real_t> cell_BM_adhesion_strength;

	ContainerType<real_t> cell_cell_repulsion_strength;
	ContainerType<real_t> cell_BM_repulsion_strength;

	ContainerType<real_t> cell_adhesion_affinities; // flattened: agents x agent_types_count

	ContainerType<real_t> relative_maximum_adhesion_distance;

	ContainerType<index_t> maximum_number_of_attachments;
	ContainerType<real_t> attachment_elastic_constant;

	ContainerType<real_t> attachment_rate;
	ContainerType<real_t> detachment_rate;

	// Motility
	ContainerType<std::uint8_t> is_motile;
	ContainerType<real_t> persistence_time;
	ContainerType<real_t> migration_speed;

	ContainerType<real_t> migration_bias_direction; // dims per agent
	ContainerType<real_t> migration_bias;

	ContainerType<real_t> motility_vector; // dims per agent

	ContainerType<std::uint8_t> restrict_to_2d;

	ContainerType<index_t> chemotaxis_index;
	ContainerType<index_t> chemotaxis_direction;
	ContainerType<real_t> chemotactic_sensitivities; // flattened: agents x substrates_count

	// State
	ContainerType<ContainerType<index_t>> neighbors;	  // neighbor indices for mechanics
	ContainerType<ContainerType<index_t>> springs;		  // spring attachments
	ContainerType<ContainerType<index_t>> attached_cells; // attachments not modeled as springs

	ContainerType<real_t> orientation;	   // dims per agent
	ContainerType<real_t> simple_pressure; // scalar mechanics pressure proxy

	ContainerType<index_t> cell_definition_indices; // link to cell definition
	ContainerType<std::uint8_t> is_movable;			// mobility toggle per agent

	explicit mechanical_agent_data(physicore::base_agent_data_generic_storage<ContainerType>& base_data,
								   index_t dims = 3, index_t agent_types_count = 0, index_t substrates_count = 0);

	void add();
	void remove_at(index_t position);
};



using agent_data = mechanical_agent_data<std::vector>;

} // namespace physicore::mechanics::physicell
