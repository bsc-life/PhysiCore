#pragma once

#include <memory>

#include <common/types.h>

namespace physicore::mechanics::micromechanics {

class environment;

/**
 * @brief Abstract interface for micromechanics solvers.
 *
 * Solvers implement the mechanical simulation algorithms for cell-cell interactions,
 * motility, membrane interactions, and position updates. Different backends (OpenMP,
 * CUDA, etc.) can provide optimized implementations.
 *
 * Following the BioFVM pattern, solvers are registered via solver_registry and
 * instantiated by name at runtime.
 */
class solver
{
public:
	/// Initialize solver with environment (allocate buffers, prepare data structures)
	virtual void initialize(environment& e) = 0;

	/// Update neighbor lists for all agents
	virtual void update_cell_neighbors(environment& e) = 0;

	/// Calculate cell-cell interaction forces using configured potentials
	virtual void update_cell_forces(environment& e) = 0;

	/// Calculate cell-level data from agents.
	/// Aggregates: positions, volumes, velocities, motility directions,
	/// compartment pressures, neighbor cells, and compartment counts.
	/// Should be called after update_cell_forces().
	virtual void calculate_cell_data(environment& e) = 0;

	/// Update motility forces based on persistence, bias, and chemotaxis
	virtual void update_motility(environment& e) = 0;

	/// Calculate basement membrane/boundary interaction forces
	virtual void update_basement_membrane_interactions(environment& e) = 0;

	/// Update spring attachment forces between connected cells
	virtual void update_spring_attachments(environment& e) = 0;

	/// Integrate positions using calculated velocities (Adams-Bashforth)
	virtual void update_positions(environment& e) = 0;

	virtual ~solver() = default;
};

using solver_ptr = std::unique_ptr<solver>;

} // namespace physicore::mechanics::micromechanics
