#pragma once

#include <cassert>
#include <cstdint>
#include <span>
#include <vector>

#include <common/base_agent_interface.h>
#include <common/types.h>

#include "cell_data.h"
#include "cell_interface.h"

namespace physicore::mechanics::micromechanics {

/**
 * @brief Template implementation of cell storage for micromechanics.
 *
 * Provides per-cell access to properties stored in SoA format in cell_data.
 */
template <typename CellDataType>
class cell_generic_storage : public virtual cell_interface
{
protected:
	CellDataType& data;

public:
	using DataType = CellDataType;
	using InterfaceType = cell_interface;

	explicit cell_generic_storage(index_t cell_id, CellDataType& data) : base_agent_interface(cell_id), data(data) {}

	// base_agent_interface
	std::span<real_t> position() override
	{
		const index_t dims = data.dims;
		assert(dims > 0);
		assert(this->index < data.cells_count);
		return std::span<real_t>(&data.positions[this->index * dims], dims);
	}

	// Geometry
	real_t& volume() override
	{
		assert(this->index < data.cells_count);
		return data.volumes[this->index];
	}

	// Definition-derived parameters / bookkeeping
	index_t& cell_definition_id() override
	{
		assert(this->index < data.cells_count);
		return data.cell_definition_ids[this->index];
	}

	index_t& compartments_count() override
	{
		assert(this->index < data.cells_count);
		return data.compartments_count[this->index];
	}

	// Kinematics
	std::span<real_t> velocity() override
	{
		const index_t dims = data.dims;
		assert(dims > 0);
		assert(this->index < data.cells_count);
		return std::span<real_t>(&data.velocities[this->index * dims], dims);
	}

	std::span<real_t> motility_direction() override
	{
		const index_t dims = data.dims;
		assert(dims > 0);
		assert(this->index < data.cells_count);
		return std::span<real_t>(&data.motility_directions[this->index * dims], dims);
	}

	// Motility configuration / state
	real_t& migration_speed() override
	{
		assert(this->index < data.cells_count);
		return data.migration_speeds[this->index];
	}

	real_t& migration_bias() override
	{
		assert(this->index < data.cells_count);
		return data.migration_biases[this->index];
	}

	std::span<real_t> polarization() override
	{
		const index_t dims = data.dims;
		assert(dims > 0);
		assert(this->index < data.cells_count);
		return std::span<real_t>(&data.polarizations[this->index * dims], dims);
	}

	std::uint8_t& compartment_is_movable(std::uint8_t agent_type) override
	{
		assert(this->index < data.cells_count);
		return data.compartment_is_movable[CellDataType::compartment_offset(this->index, agent_type)];
	}

	// Mechanics
	real_t pressure(std::uint8_t agent_type) const override
	{
		if (this->index >= data.cells_count)
			return 0.0;
		return data.compartment_pressures[CellDataType::compartment_offset(this->index, agent_type)];
	}

	void add_pressure(std::uint8_t agent_type, real_t delta) override
	{
		assert(this->index < data.cells_count);
		data.compartment_pressures[CellDataType::compartment_offset(this->index, agent_type)] += delta;
	}

	real_t total_pressure() const override
	{
		if (this->index >= data.cells_count)
			return 0.0;

		real_t total = 0.0;
		const std::size_t base = static_cast<std::size_t>(this->index) * CellDataType::compartments_per_cell;
		for (std::size_t t = 0; t < CellDataType::compartments_per_cell; ++t)
			total += data.compartment_pressures[base + t];
		return total;
	}

	index_t compartment_count(std::uint8_t agent_type) const override
	{
		if (this->index >= data.cells_count)
			return 0;
		return data.compartment_counts[CellDataType::compartment_offset(this->index, agent_type)];
	}

	index_t total_agent_count() const override
	{
		if (this->index >= data.cells_count)
			return 0;

		index_t total = 0;
		const std::size_t base = static_cast<std::size_t>(this->index) * CellDataType::compartments_per_cell;
		for (std::size_t t = 0; t < CellDataType::compartments_per_cell; ++t)
			total += data.compartment_counts[base + t];
		return total;
	}

	// Topology
	std::span<index_t> neighbor_cells() override
	{
		assert(this->index < data.cells_count);
		return std::span<index_t>(data.neighbor_cells[this->index]);
	}
};

using cell_storage = cell_generic_storage<cell_data>;

} // namespace physicore::mechanics::micromechanics
