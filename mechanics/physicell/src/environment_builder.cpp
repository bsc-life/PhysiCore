#include "environment_builder.h"

#include "solver_registry.h"
#include "vtk_agents_serializer.h"

using namespace physicore::mechanics::physicell;
using namespace physicore;

void environment_builder::set_time_step(real_t time_step) { this->timestep = time_step; }

void environment_builder::set_simulation_time(real_t sim_time) { this->simulation_time = sim_time; }

void environment_builder::resize(index_t dims, std::array<sindex_t, 3> bounding_box_mins,
								 std::array<sindex_t, 3> bounding_box_maxs, std::array<index_t, 3> voxel_shape)
{
	mesh = cartesian_mesh(dims, bounding_box_mins, bounding_box_maxs, voxel_shape);
}

void environment_builder::add_agent_type(mechanical_parameters type) { agent_types.emplace_back(std::move(type)); }
void environment_builder::add_substrate(std::string_view name) { substrate_names.emplace_back(name); }

void environment_builder::set_automated_spring_adhesion(bool value) { automated_spring_adhesion = value; }
void environment_builder::set_virtual_wall_at_domain_edges(bool value) { virtual_wall_at_domain_edges = value; }

std::unique_ptr<environment> environment_builder::build()
{
	if (!mesh)
	{
		throw std::runtime_error("Environment cannot be built without a mesh");
	}

	if (agent_types.empty())
	{
		throw std::runtime_error("Environment cannot be built with no agent types");
	}

	auto e = std::make_unique<environment>(*mesh, agent_types.size(), substrate_names.size(), timestep);

	e->agent_types = std::move(agent_types);

	e->automated_spring_adhesion = automated_spring_adhesion;
	e->virtual_wall_at_domain_edges = virtual_wall_at_domain_edges;

	e->simulation_time = simulation_time;

	auto solver = solver_registry::instance().get(solver_name);

	if (!solver)
	{
		throw std::runtime_error("Can not find solver for environment: " + solver_name);
	}

	e->solver = std::move(solver);

	e->serializer = std::make_unique<vtk_agents_serializer>("output", *e, substrate_names);

	return e;
}
