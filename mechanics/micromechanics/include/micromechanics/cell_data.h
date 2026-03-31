#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

#include <common/types.h>

namespace physicore::mechanics::micromechanics {

/**
 * @brief Cell instance data — flat SoA, one value per cell.
 *
 * In the minimal SEM framework every cell has a single compartment type
 * and all agents in a cell share the same physical properties (radius,
 * movability, motility).  Solver kernels look up these properties via
 * the agent's `cell_id`.
 *
 * Aggregate quantities (COM position/velocity, agent count, pressure)
 * are recomputed each timestep from agent data.
 */
template <template <typename...> typename ContainerType = std::vector>
struct cell_data_generic_storage
{
	static constexpr index_t invalid_cell_id = static_cast<index_t>(-1);

	index_t dims = 0;
	index_t cells_count = 0;

	// ========== Per-cell physical properties ==========
	ContainerType<real_t> radii;			// cells_count — shared agent radius
	ContainerType<std::uint8_t> is_movable; // cells_count
	ContainerType<std::uint8_t> is_motile;	// cells_count

	// ========== Motility parameters ==========
	ContainerType<real_t> migration_speeds;			 // cells_count
	ContainerType<real_t> migration_biases;			 // cells_count
	ContainerType<real_t> migration_bias_directions; // cells_count * dims
	ContainerType<real_t> persistence_times;		 // cells_count
	ContainerType<real_t> motility_directions;		 // cells_count * dims (state)

	// ========== Aggregated from agents ==========
	ContainerType<real_t> positions;	 // cells_count * dims (COM)
	ContainerType<real_t> velocities;	 // cells_count * dims (COM)
	ContainerType<index_t> agent_counts; // cells_count
	ContainerType<real_t> pressures;	 // cells_count

	// ========== Topology ==========
	ContainerType<std::vector<index_t>> neighbor_cells; // cells_count

	// ========== Helpers ==========

	/// Index into a per-cell vector field (positions, velocities, …).
	std::size_t cell_offset(index_t cell_id, index_t dim) const
	{
		return static_cast<std::size_t>(cell_id) * static_cast<std::size_t>(dims) + static_cast<std::size_t>(dim);
	}

	/// Clear all cell data vectors (count stays unchanged — call resize to reset).
	void clear();

	/// Resize / grow storage for `new_cells_count` cells in `new_dims` dimensions.
	void resize(index_t new_cells_count, index_t new_dims);

	/// Convenience overload — keeps current dims (must have been set previously).
	void resize(index_t new_cells_count);
};

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

template <template <typename...> typename ContainerType>
void cell_data_generic_storage<ContainerType>::clear()
{
	radii.clear();
	is_movable.clear();
	is_motile.clear();
	migration_speeds.clear();
	migration_biases.clear();
	migration_bias_directions.clear();
	persistence_times.clear();
	motility_directions.clear();
	positions.clear();
	velocities.clear();
	agent_counts.clear();
	pressures.clear();
	neighbor_cells.clear();
}

template <template <typename...> typename ContainerType>
void cell_data_generic_storage<ContainerType>::resize(index_t new_cells_count, index_t new_dims)
{
	assert(new_dims > 0);
	const index_t old_cells_count = cells_count;
	const index_t old_dims = dims;

	cells_count = new_cells_count;
	dims = new_dims;

	auto const n = static_cast<std::size_t>(cells_count);
	auto const nd = n * static_cast<std::size_t>(dims);

	// If dims changed, reinitialize vector-valued arrays.
	if (old_cells_count != 0 && old_dims != dims)
	{
		migration_bias_directions.clear();
		motility_directions.clear();
		positions.clear();
		velocities.clear();
	}

	// Scalar per-cell
	radii.resize(n, 1.0);
	is_movable.resize(n, static_cast<std::uint8_t>(1));
	is_motile.resize(n, static_cast<std::uint8_t>(0));
	migration_speeds.resize(n, 0.0);
	migration_biases.resize(n, 0.0);
	persistence_times.resize(n, 0.0);
	agent_counts.resize(n, 0);
	pressures.resize(n, 0.0);

	// Vector per-cell
	migration_bias_directions.resize(nd, 0.0);
	motility_directions.resize(nd, 0.0);
	positions.resize(nd, 0.0);
	velocities.resize(nd, 0.0);

	// Neighbor lists
	neighbor_cells.resize(n);
	for (index_t c = old_cells_count; c < cells_count; ++c)
		neighbor_cells[static_cast<std::size_t>(c)].clear();
}

template <template <typename...> typename ContainerType>
void cell_data_generic_storage<ContainerType>::resize(index_t new_cells_count)
{
	assert(dims > 0 && "dims must be set before calling resize(cells_count)");
	resize(new_cells_count, dims);
}

using cell_data = cell_data_generic_storage<std::vector>;

} // namespace physicore::mechanics::micromechanics
