#pragma once

#include <cstdint>
#include <span>

#include <common/base_agent_interface.h>
#include <common/types.h>

namespace physicore::mechanics::micromechanics {

/**
 * @brief Cell-level interface for micromechanics.
 *
 * Mirrors the BioFVM pattern: storage lives in cell_data, and this interface
 * exposes per-cell accessors and small operations.
 */
class cell_interface : public virtual base_agent_interface
{
public:
	// Geometry
	virtual real_t& volume() = 0;

	// Definition-derived parameters / bookkeeping
	virtual index_t& cell_definition_id() = 0;
	virtual index_t& compartments_count() = 0;

	// Kinematics
	virtual std::span<real_t> velocity() = 0;
	virtual std::span<real_t> motility_direction() = 0;

	// Motility configuration / state
	virtual real_t& migration_speed() = 0;
	virtual real_t& migration_bias() = 0;
	virtual std::span<real_t> polarization() = 0;
	virtual std::uint8_t& compartment_is_movable(std::uint8_t agent_type) = 0;

	// Mechanics (compartmented by agent_type)
	virtual real_t pressure(std::uint8_t agent_type) const = 0;
	virtual void add_pressure(std::uint8_t agent_type, real_t delta) = 0;
	virtual real_t total_pressure() const = 0;

	virtual index_t compartment_count(std::uint8_t agent_type) const = 0;
	virtual index_t total_agent_count() const = 0;

	// Topology
	virtual std::span<index_t> neighbor_cells() = 0;
};

} // namespace physicore::mechanics::micromechanics
