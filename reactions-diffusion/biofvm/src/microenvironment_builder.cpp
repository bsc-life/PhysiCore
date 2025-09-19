#include "microenvironment_builder.h"

#include <algorithm>
#include <stdexcept>

#include "microenvironment.h"
#include "types.h"

using namespace physicore::biofvm;
using namespace physicore;

void microenvironment_builder::set_name(const std::string& name) { this->name = name; }

void microenvironment_builder::set_time_units(const std::string& time_units) { this->time_units = time_units; }

void microenvironment_builder::set_space_units(const std::string& space_units) { this->space_units = space_units; }

void microenvironment_builder::set_time_step(real_t time_step) { this->timestep = time_step; }

void microenvironment_builder::resize(index_t dims, std::array<index_t, 3> bounding_box_mins,
									  std::array<index_t, 3> bounding_box_maxs, std::array<index_t, 3> voxel_shape)
{
	mesh_ = cartesian_mesh(dims, bounding_box_mins, bounding_box_maxs, voxel_shape);
}

void microenvironment_builder::add_density(const std::string& name, const std::string& units,
										   real_t diffusion_coefficient, real_t decay_rate, real_t initial_condition)
{
	substrates_names.push_back(name);
	substrates_units.push_back(units);
	diffusion_coefficients.push_back(diffusion_coefficient);
	decay_rates.push_back(decay_rate);
	initial_conditions.push_back(initial_condition);

	boundary_dirichlet_mins_values.push_back({ 0, 0, 0 });
	boundary_dirichlet_maxs_values.push_back({ 0, 0, 0 });
	boundary_dirichlet_mins_conditions.push_back({ false, false, false });
	boundary_dirichlet_maxs_conditions.push_back({ false, false, false });
}

std::size_t microenvironment_builder::get_density_index(const std::string& name) const
{
	auto it = std::find(substrates_names.begin(), substrates_names.end(), name);
	if (it == substrates_names.end())
	{
		throw std::runtime_error("Density " + name + " not found");
	}
	return std::distance(substrates_names.begin(), it);
}

void microenvironment_builder::add_dirichlet_node(std::array<index_t, 3> voxel_index, std::vector<real_t> values,
												  std::vector<bool> conditions)
{
	if (values.size() != substrates_names.size())
	{
		throw std::runtime_error("Dirichlet node values size does not match the number of densities");
	}
	if (!conditions.empty() && conditions.size() != substrates_names.size())
	{
		throw std::runtime_error("Dirichlet node conditions size does not match the number of densities");
	}
	if (!mesh_)
	{
		throw std::runtime_error("Dirichlet node cannot be added without a mesh");
	}

	if (mesh_->dims >= 1)
		dirichlet_voxels.push_back(voxel_index[0]);
	if (mesh_->dims >= 2)
		dirichlet_voxels.push_back(voxel_index[1]);
	if (mesh_->dims == 3)
		dirichlet_voxels.push_back(voxel_index[2]);

	dirichlet_values.insert(dirichlet_values.end(), values.begin(), values.end());

	if (conditions.empty())
	{
		conditions.resize(values.size(), true);
	}

	dirichlet_conditions.insert(dirichlet_conditions.end(), conditions.begin(), conditions.end());
}

void microenvironment_builder::add_boundary_dirichlet_conditions(std::size_t density_index,
																 std::array<real_t, 3> mins_values,
																 std::array<real_t, 3> maxs_values,
																 std::array<bool, 3> mins_conditions,
																 std::array<bool, 3> maxs_conditions)
{
	if (density_index >= substrates_names.size())
	{
		throw std::runtime_error("Density index out of bounds");
	}

	boundary_dirichlet_mins_values[density_index] = mins_values;
	boundary_dirichlet_maxs_values[density_index] = maxs_values;
	boundary_dirichlet_mins_conditions[density_index] = mins_conditions;
	boundary_dirichlet_maxs_conditions[density_index] = maxs_conditions;
}

void microenvironment_builder::set_bulk_functions(microenvironment::bulk_func_t supply_rate_func,
												  microenvironment::bulk_func_t uptake_rate_func,
												  microenvironment::bulk_func_t supply_target_densities_func)
{
	this->supply_rate_func = std::move(supply_rate_func);
	this->uptake_rate_func = std::move(uptake_rate_func);
	this->supply_target_densities_func = std::move(supply_target_densities_func);
}

void microenvironment_builder::do_compute_internalized_substrates() { compute_internalized_substrates = true; }

void fill_one(index_t dim_idx, index_t substrates_count, const std::vector<std::array<real_t, 3>>& values,
			  const std::vector<std::array<bool, 3>>& conditions,
			  std::array<std::unique_ptr<real_t[]>, 3>& linearized_values,
			  std::array<std::unique_ptr<bool[]>, 3>& linearized_conditions)
{
	bool any = false;
	for (index_t s = 0; s < substrates_count; s++)
	{
		any |= conditions[s][dim_idx];
	}

	if (any)
	{
		linearized_values[dim_idx] = std::make_unique<real_t[]>(substrates_count);
		linearized_conditions[dim_idx] = std::make_unique<bool[]>(substrates_count);

		for (index_t s = 0; s < substrates_count; s++)
		{
			linearized_values[dim_idx][s] = values[s][dim_idx];
			linearized_conditions[dim_idx][s] = conditions[s][dim_idx];
		}
	}
}

void microenvironment_builder::fill_dirichlet_vectors(microenvironment& m)
{
	for (index_t d = 0; d < m.mesh.dims; d++)
	{
		fill_one(d, m.substrates_count, boundary_dirichlet_mins_values, boundary_dirichlet_mins_conditions,
				 m.dirichlet_min_boundary_values, m.dirichlet_min_boundary_conditions);
		fill_one(d, m.substrates_count, boundary_dirichlet_maxs_values, boundary_dirichlet_maxs_conditions,
				 m.dirichlet_max_boundary_values, m.dirichlet_max_boundary_conditions);
	}
}

std::unique_ptr<microenvironment> microenvironment_builder::build()
{
	if (!mesh_)
	{
		throw std::runtime_error("Microenvironment cannot be built without a mesh");
	}

	if (substrates_names.empty())
	{
		throw std::runtime_error("Microenvironment cannot be built wit no densities");
	}

	auto m = std::make_unique<microenvironment>(*mesh_, substrates_names.size(), timestep);

	m->name = std::move(name);
	m->time_units = std::move(time_units);
	m->space_units = std::move(space_units);

	m->substrates_names = std::move(substrates_names);
	m->substrates_units = std::move(substrates_units);

	m->initial_conditions = std::make_unique<real_t[]>(initial_conditions.size());
	std::memcpy(m->initial_conditions.get(), initial_conditions.data(), initial_conditions.size() * sizeof(real_t));

	m->diffusion_coefficients = std::make_unique<real_t[]>(diffusion_coefficients.size());
	std::memcpy(m->diffusion_coefficients.get(), diffusion_coefficients.data(),
				diffusion_coefficients.size() * sizeof(real_t));

	m->decay_rates = std::make_unique<real_t[]>(decay_rates.size());
	std::memcpy(m->decay_rates.get(), decay_rates.data(), decay_rates.size() * sizeof(real_t));

	m->dirichlet_interior_voxels_count = dirichlet_voxels.size() / m->mesh.dims;
	m->dirichlet_interior_voxels = std::make_unique<index_t[]>(dirichlet_voxels.size());
	std::copy(dirichlet_voxels.begin(), dirichlet_voxels.end(), m->dirichlet_interior_voxels.get());

	m->dirichlet_interior_values = std::make_unique<real_t[]>(dirichlet_values.size());
	std::copy(dirichlet_values.begin(), dirichlet_values.end(), m->dirichlet_interior_values.get());

	m->dirichlet_interior_conditions = std::make_unique<bool[]>(dirichlet_conditions.size());
	std::copy(dirichlet_conditions.begin(), dirichlet_conditions.end(), m->dirichlet_interior_conditions.get());

	fill_dirichlet_vectors(*m);

	m->supply_rate_func = std::move(supply_rate_func);
	m->uptake_rate_func = std::move(uptake_rate_func);
	m->supply_target_densities_func = std::move(supply_target_densities_func);

	m->compute_internalized_substrates = compute_internalized_substrates;

	return m;
}
