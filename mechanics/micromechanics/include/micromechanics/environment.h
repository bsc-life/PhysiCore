#pragma once

#include <array>
#include <memory>

#include <common/timestep_executor.h>
#include <common/types.h>

#include "agent_container.h"
#include "cell_data.h"
#include "simulation_parameters.h"
#include "spatial_index.h"

namespace physicore::mechanics::micromechanics {

class solver;

/**
 * @brief Main environment class for micromechanics simulations.
 *
 * Holds all simulation state: agents, parameters, solver, and spatial index.
 * The solver is obtained from the solver_registry based on params.solver_name.
 */
class environment : public timestep_executor
{
public:
	/// Agent-level container (positions, velocities, forces, etc.)
	std::unique_ptr<agent_container> agents;

	/// Mechanics timestep
	real_t timestep;

	/// Cell-level data aggregated from agents (flat SoA, one value per cell)
	cell_data cells;

	/// Simulation parameters including type-based interactions
	simulation_parameters params;

	/// The solver (obtained from solver_registry)
	std::unique_ptr<solver> solver_;

	/// Spatial index for neighbor queries
	std::unique_ptr<spatial_index> index;

	/// Domain boundaries [x_min, y_min, z_min]
	std::array<real_t, 3> domain_min;

	/// Domain boundaries [x_max, y_max, z_max]
	std::array<real_t, 3> domain_max;

	explicit environment(real_t timestep);
	~environment() override;

	/**
	 * @brief Initialize the solver from the registry.
	 *
	 * Must be called after setting params.solver_name.
	 */
	void initialize_solver();

	void run_single_timestep() override;
	void serialize_state(real_t current_time) override;

	// Note: definitions live here because they are global simulation configuration
	// shared across solver backends.
};

} // namespace physicore::mechanics::micromechanics
