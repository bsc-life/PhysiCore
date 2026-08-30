#include "environment.h"

#include <memory>

#include <common/base_agent_data.h>

#include "config_reader.h"
#include "environment_builder.h"

using namespace physicore::mechanics::physicell;

environment::environment(const cartesian_mesh& mesh, index_t agent_types_count, index_t substrates_count,
						 real_t timestep)
	: mechanics_timestep(timestep), mesh(mesh), agent_types_count(agent_types_count)
{
	auto base_data = std::make_unique<physicore::base_agent_data>(mesh.dims);
	auto data = std::make_unique<mechanical_agent_data>(*base_data, agent_types_count, substrates_count);
	agents = std::make_unique<mechanical_agent_container>(std::move(base_data), std::move(data));
}

void environment::run_single_timestep() { solver->solve(*this, 1); }

void environment::serialize_state(real_t current_time)
{
	if (serializer)
	{
		serializer->serialize(*this, current_time);
	}
}

std::unique_ptr<environment> environment::create_from_config(const std::filesystem::path& config_file)
{
	// Parse the XML configuration file
	const mechanics_config config = parse_simulation_parameters(config_file);

	// Create builder
	environment_builder builder;

	// Set metadata from <overall>
	builder.set_time_step(config.overall.dt_mechanics);
	builder.set_simulation_time(config.overall.max_time);

	// Configure mesh from <domain>
	const auto& domain = config.domain;
	const index_t dims = domain.use_2D ? 2 : 3;

	const std::array<sindex_t, 3> bounding_box_mins = { static_cast<sindex_t>(domain.x_min),
														static_cast<sindex_t>(domain.y_min),
														static_cast<sindex_t>(domain.z_min) };

	const std::array<sindex_t, 3> bounding_box_maxs = { static_cast<sindex_t>(domain.x_max),
														static_cast<sindex_t>(domain.y_max),
														static_cast<sindex_t>(domain.z_max) };

	if (domain.dx <= 0 || domain.dy <= 0 || (!domain.use_2D && domain.dz <= 0))
	{
		throw std::runtime_error("Voxel dimensions must be positive");
	}

	const std::array<index_t, 3> voxel_shape = { static_cast<index_t>(domain.dx), static_cast<index_t>(domain.dy),
												 static_cast<index_t>(domain.dz) };

	builder.resize(dims, bounding_box_mins, bounding_box_maxs, voxel_shape);

	for (const auto& type : config.cell_types)
	{
		builder.add_agent_type(type.name);
	}

	// Build and return
	return builder.build();
}
