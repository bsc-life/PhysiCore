#pragma once

#include <optional>

#include <biofvm/biofvm_export.h>
#include <common/types.h>

#include "bulk_functor.h"
#include "mesh.h"
#include "microenvironment.h"

namespace physicore::reactions_diffusion::biofvm {

class BIOFVM_EXPORT microenvironment_builder
{
	std::string name, time_units, space_units;

	real_t timestep = 0.01;
	real_t simulation_time = 0.0;
	std::optional<cartesian_mesh> mesh;

	std::vector<std::string> substrates_names;
	std::vector<std::string> substrates_units;

	std::vector<real_t> diffusion_coefficients;
	std::vector<real_t> decay_rates;
	std::vector<real_t> initial_conditions;

	std::vector<index_t> dirichlet_voxels;
	std::vector<real_t> dirichlet_values;
	std::vector<bool> dirichlet_conditions;

	std::vector<std::array<real_t, 3>> boundary_dirichlet_mins_values;
	std::vector<std::array<real_t, 3>> boundary_dirichlet_maxs_values;
	std::vector<std::array<bool, 3>> boundary_dirichlet_mins_conditions;
	std::vector<std::array<bool, 3>> boundary_dirichlet_maxs_conditions;

	std::unique_ptr<bulk_functor> bulk_fnc;

	std::string solver_name = "openmp_solver";

	bool compute_internalized_substrates = false;

	void fill_dirichlet_vectors(microenvironment& m);

public:
	void set_name(std::string_view name);
	void set_time_units(std::string_view units);
	void set_space_units(std::string_view units);

	void set_time_step(real_t time_step);
	void set_simulation_time(real_t sim_time);

	// mesh functions
	void resize(index_t dims, std::array<sindex_t, 3> bounding_box_mins, std::array<sindex_t, 3> bounding_box_maxs,
				std::array<index_t, 3> voxel_shape);

	// density functions
	void add_density(const std::string& name, const std::string& units, real_t diffusion_coefficient = 0,
					 real_t decay_rate = 0, real_t initial_condition = 0);
	std::size_t get_density_index(const std::string& name) const;

	// dirichlet functions
	void add_dirichlet_node(std::array<index_t, 3> voxel_index, std::vector<real_t> values,
							std::vector<bool> conditions = {});
	void add_boundary_dirichlet_conditions(std::size_t density_index, std::array<real_t, 3> mins_values,
										   std::array<real_t, 3> maxs_values,
										   std::array<bool, 3> mins_conditions = { true, true, true },
										   std::array<bool, 3> maxs_conditions = { true, true, true });

	void set_bulk_functions(std::unique_ptr<bulk_functor> bulk_fnc);

	void do_compute_internalized_substrates();

	void select_solver(const std::string& solver_name);

	std::unique_ptr<microenvironment> build();
};

} // namespace physicore::reactions_diffusion::biofvm
