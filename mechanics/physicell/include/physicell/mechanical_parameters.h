#pragma once

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

enum class chemotaxis_direction_kind
{
	TOWARD_GRADIENT,
	AWAY_FROM_GRADIENT,
	NONE
};

struct mechanical_parameters
{
	// Cell type identification
	index_t id;
	std::string name;

	// Geometry parameters
	real_t radius;

	// Mechanics parameters
	real_t cell_cell_adhesion_strength;
	real_t cell_cell_repulsion_strength;
	real_t relative_maximum_adhesion_distance;
	std::vector<real_t> cell_adhesion_affinity;
	index_t maximum_number_of_attachments;
	real_t cell_BM_adhesion_strength;
	real_t cell_BM_repulsion_strength;
	real_t attachment_elastic_coefficient;
	real_t attachment_rate;
	real_t detachment_rate;

	// Motility parameters
	real_t motility_speed;
	real_t motility_persistence_time;
	real_t motility_bias;
	bool is_motile;
	bool use_2D;

	// Chemotaxis parameters
	chemotaxis_direction_kind simple_chemotaxis_direction;
	index_t simple_chemotaxis_substrate;
	bool advanced_chemotaxis_enabled;
	bool advanced_chemotaxis_normalize_each_gradient;
	std::vector<real_t> chemotaxis_sensitivities;
};


} // namespace physicore::mechanics::physicell
