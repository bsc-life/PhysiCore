#include "micromechanics/environment.h"

#include "micromechanics/solver.h"
#include "micromechanics/solver_registry.h"
#include "micromechanics/uniform_grid_spatial_index.h"
#include "micromechanics/vtk_mechanics_serializer.h"

namespace physicore::mechanics::micromechanics {

environment::environment(real_t timestep) : timestep(timestep)
{
	auto base_data = std::make_unique<physicore::base_agent_data>(3);
	auto mech_data = std::make_unique<agent_data>(*base_data);
	agents = std::make_unique<agent_container>(std::move(base_data), std::move(mech_data));
	index = std::make_unique<uniform_grid_spatial_index>();
	serializer = std::make_unique<vtk_mechanics_serializer>("output");
}

environment::~environment() = default;

void environment::initialize_solver()
{
	active_solver = solver_registry::instance().get(params.solver_name);
	if (active_solver)
	{
		active_solver->initialize(*this);
	}
}

void environment::run_single_timestep()
{
	if (active_solver)
	{
		active_solver->update_cell_neighbors(*this);
		active_solver->update_cell_forces(*this);

		// Calculate cell-level data (positions, volumes, velocities, pressures, neighbors, etc.)
		active_solver->calculate_cell_data(*this);

		if (params.enable_motility)
		{
			active_solver->update_motility(*this);
		}

		if (params.enable_basement_membrane)
		{
			active_solver->update_basement_membrane_interactions(*this);
		}

		if (params.enable_spring_attachments)
		{
			active_solver->update_spring_attachments(*this);
		}

		active_solver->update_positions(*this);
	}
}

void environment::serialize_state(real_t current_time)
{
	if (serializer)
	{
		serializer->serialize(*this, current_time);
	}
}

} // namespace physicore::mechanics::micromechanics
