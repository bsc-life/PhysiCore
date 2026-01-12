#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <common/types.h>

/*
	mechanical_parameters
	Mechanics tunables for a single cell type, mirroring the common PhysiCell "mechanics" XML entries.
	Contains the cell type ID, name, and all mechanical parameters for that cell type.
	These are intended to be stored in a container (e.g. vector<mechanical_parameters>)
	and kept distinct from global timestep / executor settings (dt, agent_types_count).
*/
namespace physicore::mechanics::physicell {

struct mechanical_parameters
{
	// Cell type identification
	index_t id;
	std::string name;

	// Mechanics parameters
	real_t cell_cell_adhesion_strength;
	real_t cell_cell_repulsion_strength;
	real_t relative_maximum_adhesion_distance;

	// Adhesion affinity to other cell types (indexed by other cell type ID)
	std::vector<real_t> cell_adhesion_affinity;

	// Options
	real_t set_relative_maximum_adhesion_distance;
	real_t set_absolute_maximum_adhesion_distance;

	index_t maximum_number_of_attachments;

	real_t cell_BM_adhesion_strength;
	real_t cell_BM_repulsion_strength;

	// Attachment parameters
	real_t attachment_elastic_coefficient;
	real_t attachment_rate;
	real_t detachment_rate;

	// Motility parameters
	real_t motility_speed;
	real_t motility_persistence_time;
	real_t motility_bias;
	bool is_movable;

	// Chemotaxis parameters (indexed by substrate ID)
	std::vector<real_t> chemotaxis_sensitivity;
	std::vector<bool> chemotaxis_enabled;

	// Advanced motility
	bool normalize_each_gradient;
	std::vector<real_t> chemotaxis_advanced_enabled;
};


} // namespace physicore::mechanics::physicell
