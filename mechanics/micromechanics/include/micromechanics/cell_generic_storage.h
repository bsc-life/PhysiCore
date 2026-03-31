#pragma once

#include <cassert>
#include <span>

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
		return std::span<real_t>(&data.positions[data.cell_offset(this->index, 0)], static_cast<std::size_t>(dims));
	}

	// Per-cell physical property
	real_t& radius() override
	{
		assert(this->index < data.cells_count);
		return data.radii[static_cast<std::size_t>(this->index)];
	}

	// Kinematics
	std::span<real_t> velocity() override
	{
		const index_t dims = data.dims;
		assert(dims > 0);
		assert(this->index < data.cells_count);
		return std::span<real_t>(&data.velocities[data.cell_offset(this->index, 0)], static_cast<std::size_t>(dims));
	}

	// Motility
	real_t& migration_speed() override
	{
		assert(this->index < data.cells_count);
		return data.migration_speeds[static_cast<std::size_t>(this->index)];
	}

	real_t& migration_bias() override
	{
		assert(this->index < data.cells_count);
		return data.migration_biases[static_cast<std::size_t>(this->index)];
	}

	// Aggregated quantities
	real_t pressure() const override
	{
		if (this->index >= data.cells_count)
			return 0.0;
		return data.pressures[static_cast<std::size_t>(this->index)];
	}

	index_t agent_count() const override
	{
		if (this->index >= data.cells_count)
			return 0;
		return data.agent_counts[static_cast<std::size_t>(this->index)];
	}

	// Topology
	std::span<index_t> neighbor_cells() override
	{
		assert(this->index < data.cells_count);
		return std::span<index_t>(data.neighbor_cells[static_cast<std::size_t>(this->index)]);
	}
};

using cell_storage = cell_generic_storage<cell_data>;

} // namespace physicore::mechanics::micromechanics
