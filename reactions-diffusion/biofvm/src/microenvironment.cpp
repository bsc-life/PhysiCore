#include "microenvironment.h"

#include "agent_container.h"
#include "base_agent_data.h"
#include "config_reader.h"
#include "microenvironment_builder.h"

using namespace physicore;
using namespace physicore::biofvm;

microenvironment::microenvironment(const cartesian_mesh& mesh, index_t substrates_count, real_t timestep)
	: diffusion_timestep(timestep), mesh(mesh), substrates_count(substrates_count)
{
	auto base_data = std::make_unique<base_agent_data>(mesh.dims);
	auto data = std::make_unique<agent_data>(*base_data, substrates_count);
	agents = make_unique<agent_container>(std::move(base_data), std::move(data));
}

std::unique_ptr<microenvironment> microenvironment::create_from_config(const std::filesystem::path& config_file)
{
	// Parse the XML configuration file
	physicell_config config = parse_physicell_config(config_file);

	// Create builder
	microenvironment_builder builder;

	// Set metadata from <overall>
	builder.set_name("microenvironment");
	builder.set_time_units(config.overall.time_units);
	builder.set_space_units(config.overall.space_units);
	builder.set_time_step(config.overall.dt_diffusion);
	builder.set_simulation_time(config.overall.max_time);

	// Configure mesh from <domain>
	const auto& domain = config.domain;
	index_t dims = domain.use_2D ? 2 : 3;

	std::array<sindex_t, 3> bounding_box_mins = { static_cast<sindex_t>(domain.x_min),
												  static_cast<sindex_t>(domain.y_min),
												  static_cast<sindex_t>(domain.z_min) };

	std::array<sindex_t, 3> bounding_box_maxs = { static_cast<sindex_t>(domain.x_max),
												  static_cast<sindex_t>(domain.y_max),
												  static_cast<sindex_t>(domain.z_max) };

	if (domain.dx <= 0 || domain.dy <= 0 || (!domain.use_2D && domain.dz <= 0))
	{
		throw std::runtime_error("Voxel dimensions must be positive");
	}

	std::array<index_t, 3> voxel_shape = { static_cast<index_t>(domain.dx), static_cast<index_t>(domain.dy),
										   static_cast<index_t>(domain.dz) };

	builder.resize(dims, bounding_box_mins, bounding_box_maxs, voxel_shape);

	// Add substrates from <microenvironment_setup>
	for (const auto& variable : config.microenvironment.variables)
	{
		builder.add_density(variable.name, variable.units, variable.diffusion_coefficient, variable.decay_rate,
							variable.initial_condition);

		// Add boundary Dirichlet conditions for this substrate
		std::size_t density_index = builder.get_density_index(variable.name);
		builder.add_boundary_dirichlet_conditions(
			density_index, variable.boundary_conditions.mins_values, variable.boundary_conditions.maxs_values,
			variable.boundary_conditions.mins_conditions, variable.boundary_conditions.maxs_conditions);
	}

	// Set options
	if (config.microenvironment.track_internalized_substrates)
	{
		builder.do_compute_internalized_substrates();
	}

	// Build and return
	return builder.build();
}

void microenvironment::run_single_timestep() { solver->solve(*this, 1); }

void microenvironment::serialize_state(real_t current_time)
{
	if (serializer)
		serializer->serialize(*this, current_time);
	if (agents_serializer)
		agents_serializer->serialize(*this, current_time);
}

real_t microenvironment::get_substrate_density(index_t s, index_t x, index_t y, index_t z) const
{
	return solver->get_substrate_density(s, x, y, z);
}

void microenvironment::print_info(std::ostream& os) const
{
	os << "Microenvironment config:" << std::endl;
	os << "  Time units: " << time_units << std::endl;
	os << "  Space units: " << space_units << std::endl;
	os << "  Timestep: " << diffusion_timestep << " " << time_units << std::endl;
	os << "  Dimensions: " << mesh.dims << "D" << std::endl;
	os << "  Grid bounds: [" << mesh.bounding_box_mins[0] << ", " << mesh.bounding_box_maxs[0] << "] x ["
	   << mesh.bounding_box_mins[1] << ", " << mesh.bounding_box_maxs[1] << "]";
	if (mesh.dims == 3)
	{
		os << " x [" << mesh.bounding_box_mins[2] << ", " << mesh.bounding_box_maxs[2] << "]";
	}
	os << " " << space_units << std::endl;
	os << "  Voxel size: " << mesh.voxel_shape[0] << " x " << mesh.voxel_shape[1];
	if (mesh.dims == 3)
	{
		os << " x " << mesh.voxel_shape[2];
	}
	os << " " << space_units << std::endl;
	os << "  Grid resolution: " << mesh.grid_shape[0] << " x " << mesh.grid_shape[1];
	if (mesh.dims == 3)
	{
		os << " x " << mesh.grid_shape[2];
	}
	os << " voxels" << std::endl;
	os << "  Substrates: " << substrates_count << std::endl;

	for (index_t i = 0; i < substrates_count; ++i)
	{
		os << "    - " << substrates_names[i] << " (D=" << diffusion_coefficients[i] << ", λ=" << decay_rates[i]
		   << ", I=" << initial_conditions[i] << ")" << std::endl;
	}
}
