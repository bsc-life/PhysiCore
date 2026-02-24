#pragma once

#include <common/types.h>

namespace physicore::mechanics::physicell::kernels::openmp_solver {

/**
 * @brief Pairwise cell-cell force calculation interface
 *
 * Computes repulsive and adhesive forces between two cells based on:
 * - Overlap (repulsive potential)
 * - Affinity-based adhesion (within adhesion distance)
 * - Cell radii and mechanical properties
 *
 * Implements Newton's third law: force on cell lhs from rhs is equal
 * and opposite to force on cell rhs from lhs.
 *
 * @tparam dims Dimensionality (1, 2, or 3)
 * @param lhs Index of first cell
 * @param rhs Index of second cell
 * @param cell_defs_count Number of cell types
 * @param velocity Velocity array (dims * num_cells)
 * @param simple_pressure Simple pressure array (accumulated repulsion metric)
 * @param position Position array (dims * num_cells)
 * @param radius Cell radius array (num_cells)
 * @param cell_cell_repulsion_strength Repulsion strength per cell (num_cells)
 * @param cell_cell_adhesion_strength Adhesion strength per cell (num_cells)
 * @param relative_maximum_adhesion_distance Adhesion range multiplier per cell (num_cells)
 * @param cell_adhesion_affinity Affinity matrix for cell types (cell_defs_count * cell_defs_count)
 * @param cell_definition_index Cell type index per cell (num_cells)
 */
template <physicore::index_t dims>
void solve_pair(physicore::index_t lhs, physicore::index_t rhs, physicore::index_t cell_defs_count,
				physicore::real_t* __restrict__ velocity, physicore::real_t* __restrict__ simple_pressure,
				const physicore::real_t* __restrict__ position, const physicore::real_t* __restrict__ radius,
				const physicore::real_t* __restrict__ cell_cell_repulsion_strength,
				const physicore::real_t* __restrict__ cell_cell_adhesion_strength,
				const physicore::real_t* __restrict__ relative_maximum_adhesion_distance,
				const physicore::real_t* __restrict__ cell_adhesion_affinity,
				const physicore::index_t* __restrict__ cell_definition_index);

} // namespace physicore::mechanics::physicell::kernels::openmp_solver
