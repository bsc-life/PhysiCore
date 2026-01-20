#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

#include <common/types.h>

namespace physicore::mechanics::micromechanics {

/**
 * @brief Cell instance data - runtime state for each cell.
 *
 * All properties are indexed by cell_id. Properties are computed
 * by aggregating data from agents belonging to each cell.
 */
template <template <typename...> typename ContainerType = std::vector>
struct cell_data_generic_storage
{
	// Matches agent_data design: flat SoA vectors + explicit indexing.
	// Assumption (as per reviewer request): cell_id is dense enough to index directly.
	static constexpr index_t invalid_cell_id = static_cast<index_t>(-1);
	static constexpr index_t default_dims = 3;
	static constexpr index_t compartments_per_cell = 256; // agent_type is uint8_t

	index_t cells_count = 0;
	index_t dims = default_dims;

	// ========== Geometry ==========
	ContainerType<real_t> positions; // cells_count * dims
	ContainerType<real_t> volumes;	 // cells_count

	// ========== Cell definition bookkeeping ==========
	// Index into a (future) cell-definition table.
	ContainerType<index_t> cell_definition_ids; // cells_count
	// Number of compartments for this cell (derived from its definition).
	ContainerType<index_t> compartments_count; // cells_count

	// ========== Kinematics ==========
	ContainerType<real_t> velocities;		   // cells_count * dims
	ContainerType<real_t> motility_directions; // cells_count * dims

	// ========== Motility configuration / state ==========
	// PhysiCell-like parameters, defined per cell (often derived from a cell definition).
	ContainerType<real_t> migration_speeds; // cells_count
	ContainerType<real_t> migration_biases; // cells_count
	// Persistent random walk / polarization vector (unit or not, model-dependent).
	ContainerType<real_t> polarizations; // cells_count * dims
	// Per-cell per-compartment movable flag (0/1), indexed by agent_type.
	ContainerType<std::uint8_t> compartment_is_movable; // cells_count * 256

	// ========== Mechanics ==========
	ContainerType<real_t> compartment_pressures; // cells_count * 256
	ContainerType<index_t> compartment_counts;	 // cells_count * 256

	// ========== Topology ==========
	ContainerType<std::vector<index_t>> neighbor_cells; // cells_count

	// ========== Methods ==========

	/// Clear all cell data.
	void clear();

	/// Resize storage for a given number of cells (zero-initialize new entries).
	void resize(index_t new_cells_count, index_t new_dims = default_dims);

	constexpr std::size_t cell_offset(index_t cell_id, index_t dim) const
	{
		return static_cast<std::size_t>(cell_id) * static_cast<std::size_t>(dims) + static_cast<std::size_t>(dim);
	}
	static constexpr std::size_t compartment_offset(index_t cell_id, std::uint8_t agent_type)
	{
		return static_cast<std::size_t>(cell_id) * static_cast<std::size_t>(compartments_per_cell)
			   + static_cast<std::size_t>(agent_type);
	}
};

template <template <typename...> typename ContainerType>
void cell_data_generic_storage<ContainerType>::clear()
{
	positions.clear();
	volumes.clear();
	cell_definition_ids.clear();
	compartments_count.clear();
	velocities.clear();
	motility_directions.clear();
	migration_speeds.clear();
	migration_biases.clear();
	polarizations.clear();
	compartment_is_movable.clear();
	compartment_pressures.clear();
	compartment_counts.clear();
	neighbor_cells.clear();

	cells_count = 0;
	dims = default_dims;
}

template <template <typename...> typename ContainerType>
void cell_data_generic_storage<ContainerType>::resize(index_t new_cells_count, index_t new_dims)
{
	assert(new_dims > 0);
	const index_t old_cells_count = cells_count;
	const index_t old_dims = dims;

	cells_count = new_cells_count;
	dims = new_dims;

	// If dims changed, simplest is to reinitialize geometry/kinematics arrays.
	if (old_cells_count != 0 && old_dims != dims)
	{
		positions.clear();
		velocities.clear();
		motility_directions.clear();
		polarizations.clear();
	}

	positions.resize(static_cast<std::size_t>(cells_count) * static_cast<std::size_t>(dims), 0.0);
	velocities.resize(static_cast<std::size_t>(cells_count) * static_cast<std::size_t>(dims), 0.0);
	motility_directions.resize(static_cast<std::size_t>(cells_count) * static_cast<std::size_t>(dims), 0.0);
	polarizations.resize(static_cast<std::size_t>(cells_count) * static_cast<std::size_t>(dims), 0.0);

	volumes.resize(static_cast<std::size_t>(cells_count), 0.0);
	cell_definition_ids.resize(static_cast<std::size_t>(cells_count), 0);
	compartments_count.resize(static_cast<std::size_t>(cells_count), 0);
	migration_speeds.resize(static_cast<std::size_t>(cells_count), 0.0);
	migration_biases.resize(static_cast<std::size_t>(cells_count), 0.0);

	compartment_pressures.resize(static_cast<std::size_t>(cells_count) * compartments_per_cell, 0.0);
	compartment_counts.resize(static_cast<std::size_t>(cells_count) * compartments_per_cell, 0);
	compartment_is_movable.resize(static_cast<std::size_t>(cells_count) * compartments_per_cell,
								  static_cast<std::uint8_t>(1));

	neighbor_cells.resize(static_cast<std::size_t>(cells_count));
	// On growth, clear only the newly created neighbor vectors.
	for (index_t c = old_cells_count; c < cells_count; ++c)
		neighbor_cells[static_cast<std::size_t>(c)].clear();
}

using cell_data = cell_data_generic_storage<std::vector>;

} // namespace physicore::mechanics::micromechanics
