#pragma once
#include <cstdint>
#include <span>

#include <common/base_agent_interface.h>
#include <common/types.h>

namespace physicore::mechanics::micromechanics {

/**
 * @brief Interface for mechanics agents (sub-cellular compartments).
 *
 * This interface extends base_agent_interface with mechanics-specific
 * properties. Note that cell-level properties (pressure, volume, etc.)
 * are accessed via cell_data, not through this interface.
 *
 * Potential parameters and interaction configs are accessed via
 * simulation_parameters at the cell/type level.
 */
class cell_interface : public virtual base_agent_interface
{
public:
	// Kinematics
	virtual std::span<real_t> velocity() = 0;
	virtual std::span<real_t> previous_velocity() = 0;

	// Geometry
	virtual real_t& radius() = 0;
	virtual std::uint8_t& is_movable() = 0;

	// Per-agent interaction strengths (can vary within a cell type)
	virtual real_t& cell_cell_adhesion_strength() = 0;
	virtual real_t& cell_cell_repulsion_strength() = 0;
	virtual real_t& relative_maximum_adhesion_distance() = 0;

	// Topology
	virtual std::span<index_t> neighbors() = 0;
};

} // namespace physicore::mechanics::micromechanics
