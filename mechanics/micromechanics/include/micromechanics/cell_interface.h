#pragma once

#include <span>

#include <common/base_agent_interface.h>
#include <common/types.h>

namespace physicore::mechanics::micromechanics {

/**
 * @brief Cell-level interface for micromechanics (minimal SEM).
 *
 * Storage lives in cell_data; this interface exposes per-cell accessors.
 */
class cell_interface : public virtual base_agent_interface
{
public:
	// Per-cell physical property
	virtual real_t& radius() = 0;

	// Kinematics
	virtual std::span<real_t> velocity() = 0;

	// Motility
	virtual real_t& migration_speed() = 0;
	virtual real_t& migration_bias() = 0;

	// Aggregated quantities
	virtual real_t pressure() const = 0;
	virtual index_t agent_count() const = 0;

	// Topology
	virtual std::span<index_t> neighbor_cells() = 0;
};

} // namespace physicore::mechanics::micromechanics
