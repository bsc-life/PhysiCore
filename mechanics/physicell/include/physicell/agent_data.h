#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <common/base_agent_data.h>
#include <common/types.h>

namespace physicore::mechanics::physicell {

struct mechanical_parameters;

/**
 * @brief Mechanics-related properties for agents
 * Includes adhesion, repulsion, and attachment mechanics
 */
struct mechanics_properties
{
	// Adhesion forces
	std::vector<real_t> cell_cell_adhesion_strength;
	std::vector<real_t> cell_BM_adhesion_strength;

	// Repulsion forces
	std::vector<real_t> cell_cell_repulsion_strength;
	std::vector<real_t> cell_BM_repulsion_strength;

	// Affinity-based interaction (per agent type)
	std::vector<real_t> cell_adhesion_affinities; // flattened: agents x agent_types_count
	std::vector<real_t> relative_maximum_adhesion_distance;

	// Attachment mechanics
	std::vector<index_t> maximum_number_of_attachments;
	std::vector<real_t> attachment_elastic_constant;
	std::vector<real_t> attachment_rate;
	std::vector<real_t> detachment_rate;
};

/**
 * @brief Motility-related properties for agents
 * Includes migration speed, persistence, and chemotaxis
 */
struct motility_properties
{
	// Migration parameters
	std::vector<std::uint8_t> is_motile;
	std::vector<real_t> persistence_time;
	std::vector<real_t> migration_speed;
	std::vector<real_t> migration_bias_direction; // dims per agent
	std::vector<real_t> migration_bias;
	std::vector<real_t> motility_vector; // dims per agent
	std::vector<std::uint8_t> restrict_to_2d;

	// Chemotaxis parameters
	std::vector<index_t> chemotaxis_index;
	std::vector<index_t> chemotaxis_direction;
	std::vector<real_t> chemotactic_sensitivities; // flattened: agents x substrates_count
};

/**
 * @brief State-related properties for agents
 * Includes neighbor tracking, attachments, orientation, and mobility
 */
struct state_properties
{
	// Spatial relationships
	std::vector<std::vector<index_t>> neighbors;	  // neighbor indices for mechanics
	std::vector<std::vector<index_t>> springs;		  // spring attachments
	std::vector<std::vector<index_t>> attached_cells; // attachments not modeled as springs

	// Orientation and pressure
	std::vector<real_t> orientation;	   // dims per agent
	std::vector<real_t> simple_pressure; // scalar mechanics pressure proxy

	// Cell metadata
	std::vector<index_t> agent_type_index; // link to cell definition
	std::vector<std::uint8_t> is_movable;			// mobility toggle per agent
};

struct mechanical_agent_data
{
public:
	physicore::base_agent_data_generic_storage<std::vector>& base_data;

	// Metadata
	index_t agents_count = 0;
	index_t agent_types_count = 1;
	index_t substrates_count = 1;

	// Geometry / kinematics (separate from mechanics/motility/state)
	std::vector<real_t> velocity;		   // current velocity per axis
	std::vector<real_t> previous_velocity; // previous velocity per axis
	std::vector<real_t> radius;			   // cell radius

	// Organized sub-structures
	mechanics_properties mechanics;
	motility_properties motility;
	state_properties state;

	explicit mechanical_agent_data(physicore::base_agent_data_generic_storage<std::vector>& base_data,
								   index_t agent_types_count = 1, index_t substrates_count = 1);

	void add();
	void remove_at(index_t position);
	void add(index_t new_id, index_t cell_type,mechanical_parameters& parameters, bool is_2D);
};

using agent_data = mechanical_agent_data;

} // namespace physicore::mechanics::physicell
