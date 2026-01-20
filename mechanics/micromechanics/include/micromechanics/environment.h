#pragma once

#include <array>
#include <memory>

#include <common/timestep_executor.h>
#include <common/types.h>

#include "cell_data.h"
#include "simulation_parameters.h"
#include "spatial_index.h"

namespace physicore::mechanics::micromechanics {

class agent_container;
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
	/// Mechanics timestep
	real_t timestep;

	/// Agent container with all agent data
	std::unique_ptr<agent_container> agents;

	/// Cell-level data (pressure, etc.) aggregated from agents
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

	// TODO cell definition class + vector of cell definitions
};

} // namespace physicore::mechanics::micromechanics
