#include "micromechanics/environment.h"

#include "micromechanics/agent_container.h"
#include "micromechanics/solver.h"
#include "micromechanics/solver_registry.h"
#include "micromechanics/uniform_grid_spatial_index.h"

namespace physicore::mechanics::micromechanics {

environment::environment(real_t timestep) : timestep(timestep)
{
	auto base_data = std::make_unique<physicore::base_agent_data>(3);
	auto mech_data = std::make_unique<agent_data>(*base_data);
	agents = std::make_unique<agent_container>(std::move(base_data), std::move(mech_data));
	index = std::make_unique<uniform_grid_spatial_index>();

	domain_min = { -500.0, -500.0, -500.0 };
	domain_max = { 500.0, 500.0, 500.0 };
}

environment::~environment() = default;

void environment::initialize_solver()
{
	solver_ = solver_registry::instance().get(params.solver_name);
	if (solver_)
	{
		solver_->initialize(*this);
	}
}

void environment::run_single_timestep()
{
	if (solver_)
	{
		solver_->update_cell_neighbors(*this);
		solver_->update_cell_forces(*this);

		// Calculate cell-level data (positions, volumes, velocities, pressures, neighbors, etc.)
		solver_->calculate_cell_data(*this);

		if (params.enable_motility)
		{
			solver_->update_motility(*this);
		}

		if (params.enable_basement_membrane)
		{
			solver_->update_basement_membrane_interactions(*this);
		}

		if (params.enable_spring_attachments)
		{
			solver_->update_spring_attachments(*this);
		}

		solver_->update_positions(*this);
	}
}

void environment::serialize_state(real_t current_time) { (void)current_time; }

} // namespace physicore::mechanics::micromechanics
