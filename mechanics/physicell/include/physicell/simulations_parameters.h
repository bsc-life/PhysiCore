// File: mechanics/physicell/include/simulations_parameters.h
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <common/types.h>

/*
	SimulationParameters
	Per-agent-type mechanics tunables mirroring the common PhysiCell "mechanics" XML entries.
	These are intended to be stored in environment (e.g. environment.h: vector<SimulationParameters>)
	and kept distinct from global timestep / executor settings (dt, agent_types_count).
*/
namespace physicore::mechanics::physicell {

template <template <typename...> typename ContainerType = std::vector>
struct SimulationParameters
{
	ContainerType<bool> is_movable; // Per cell type

	// From xml
	//  Mechanics parameters
	//  Adhesion parameters
	ContainerType<real_t> cell_cell_adhesion_strength;	// Per cell type
	ContainerType<real_t> cell_cell_repulsion_strength; // Per cell type

	// Adhesion parameters
	ContainerType<real_t> relative_maximum_adhesion_distance;	 // Per cell type
	ContainerType<ContainerType<real_t>> cell_adhesion_affinity; // Cell type x Cell type

	// Options can be stored in a map if needed
	ContainerType<real_t> set_relative_maximum_adhesion_distance; // Per cell type
	ContainerType<real_t> set_absolute_maximum_adhesion_distance; // Per cell type

	ContainerType<index_t> maximum_number_of_attachments; // Per cell type

	ContainerType<real_t> cell_BM_adhesion_strength;  // Per cell type
	ContainerType<real_t> cell_BM_repulsion_strength; // Per cell type

	ContainerType<real_t> attachment_elastic_coefficient; // Per cell type
	ContainerType<real_t> attachment_rate;				  // Per cell type
	ContainerType<real_t> detachment_rate;				  // Per cell type

	// Motility parameters
	ContainerType<real_t> motility_speed;			 // Per cell type
	ContainerType<real_t> motility_persistence_time; // Per cell type
	ContainerType<real_t> motility_bias;			 // Per cell type

	// Basic motility
	bool basic_motility_enabled;
	bool is_2D;
	ContainerType<ContainerType<real_t>> chemotaxis_sensitivity; // Cell type x Substrate type
	ContainerType<ContainerType<bool>> chemotaxis_enabled;		 // Cell type x Substrate type

	// Advanced motility -
	bool advanced_motility_enabled;
	ContainerType<bool> normalize_each_gradient;					  // Per cell type
	ContainerType<ContainerType<real_t>> chemotaxis_advanced_enabled; // Cell type x Substrate type
};

} // namespace physicore::mechanics::physicell
