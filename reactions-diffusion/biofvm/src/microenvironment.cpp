#include "microenvironment.h"

#include <algorithm>

#include <common/base_agent_data.h>

#include "agent_container.h"
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

	if (!config.solver.name.empty())
	{
		builder.select_solver(config.solver.name);
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

real_t& microenvironment::get_substrate_density(index_t s, index_t x, index_t y, index_t z)
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

void microenvironment::update_dirichlet_interior_voxel(std::array<index_t, 3> voxel, index_t substrate_idx,
													   real_t value, bool condition)
{
	if (substrate_idx >= substrates_count)
	{
		throw std::runtime_error("Substrate index out of bounds");
	}

	// Try to find the voxel in the dirichlet_interior_voxels array
	for (index_t i = 0; i < dirichlet_interior_voxels_count; ++i)
	{
		auto interior_voxel_idx = &dirichlet_interior_voxels[i * mesh.dims];
		if (mesh.dims >= 1 && voxel[0] != interior_voxel_idx[0])
			continue;
		if (mesh.dims >= 2 && voxel[1] != interior_voxel_idx[1])
			continue;
		if (mesh.dims >= 3 && voxel[2] != interior_voxel_idx[2])
			continue;

		// Found the voxel, update value and condition
		index_t offset = i * substrates_count + substrate_idx;
		dirichlet_interior_values[offset] = value;
		dirichlet_interior_conditions[offset] = condition;

		return;
	}

	// Voxel not found, need to add a new entry
	index_t new_count = dirichlet_interior_voxels_count + 1;
	auto new_voxels = std::make_unique<index_t[]>(new_count * mesh.dims);
	auto new_values = std::make_unique<real_t[]>(new_count * substrates_count);
	auto new_conditions = std::make_unique<bool[]>(new_count * substrates_count);

	// Copy old data
	std::copy(dirichlet_interior_voxels.get(),
			  dirichlet_interior_voxels.get() + dirichlet_interior_voxels_count * mesh.dims, new_voxels.get());
	std::copy(dirichlet_interior_values.get(),
			  dirichlet_interior_values.get() + dirichlet_interior_voxels_count * substrates_count, new_values.get());
	std::copy(dirichlet_interior_conditions.get(),
			  dirichlet_interior_conditions.get() + dirichlet_interior_voxels_count * substrates_count,
			  new_conditions.get());

	// Add new voxel
	std::copy(voxel.begin(), voxel.begin() + mesh.dims, &new_voxels[dirichlet_interior_voxels_count * mesh.dims]);

	// Initialize new values and conditions
	for (index_t s = 0; s < substrates_count; ++s)
	{
		index_t offset = dirichlet_interior_voxels_count * substrates_count + s;
		if (s == substrate_idx)
		{
			new_values[offset] = value;
			new_conditions[offset] = condition;
		}
		else
		{
			new_values[offset] = 0.0;
			new_conditions[offset] = false;
		}
	}

	// Update pointers and count
	dirichlet_interior_voxels = std::move(new_voxels);
	dirichlet_interior_values = std::move(new_values);
	dirichlet_interior_conditions = std::move(new_conditions);
	dirichlet_interior_voxels_count = new_count;
}

namespace {
void update_boundary(microenvironment& m, char dimension, index_t substrate_idx, real_t value, bool condition,
					 auto& boundary_values, auto& boundary_conditions)
{
	index_t dim_idx = static_cast<index_t>(dimension - 'x');

	if (dim_idx >= m.mesh.dims)
	{
		throw std::runtime_error("Dimension index out of bounds");
	}

	if (substrate_idx >= m.substrates_count)
	{
		throw std::runtime_error("Substrate index out of bounds");
	}

	if (!boundary_values[dim_idx])
	{
		// Allocate arrays
		boundary_values[dim_idx] = std::make_unique<real_t[]>(m.substrates_count);
		boundary_conditions[dim_idx] = std::make_unique<bool[]>(m.substrates_count);

		// Initialize to defaults
		for (index_t s = 0; s < m.substrates_count; ++s)
		{
			boundary_values[dim_idx][s] = 0.0;
			boundary_conditions[dim_idx][s] = false;
		}
	}

	boundary_values[dim_idx][substrate_idx] = value;
	boundary_conditions[dim_idx][substrate_idx] = condition;
}
} // namespace

void microenvironment::update_dirichlet_boundary_min(char dimension, index_t substrate_idx, real_t value,
													 bool condition)
{
	update_boundary(*this, dimension, substrate_idx, value, condition, dirichlet_min_boundary_values,
					dirichlet_min_boundary_conditions);
}

void microenvironment::update_dirichlet_boundary_max(char dimension, index_t substrate_idx, real_t value,
													 bool condition)
{
	update_boundary(*this, dimension, substrate_idx, value, condition, dirichlet_max_boundary_values,
					dirichlet_max_boundary_conditions);
}

void microenvironment::update_dirichlet_conditions() { solver->reinitialize_dirichlet(*this); }
