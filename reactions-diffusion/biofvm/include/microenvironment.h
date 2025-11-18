#pragma once

#include <filesystem>
#include <memory>
#include <string>
#include <timestep_executor.h>
#include <types.h>
#include <vector>

#include "agent_container.h"
#include "bulk_functor.h"
#include "mesh.h"
#include "serializer.h"
#include "solver.h"

namespace physicore::biofvm {

class microenvironment : public timestep_executor
{
public:
	microenvironment(const cartesian_mesh& mesh, index_t substrates_count, real_t timestep);

	microenvironment(microenvironment&&) = delete;
	microenvironment(const microenvironment&) = delete;
	microenvironment& operator=(const microenvironment&) = delete;
	microenvironment& operator=(microenvironment&&) = delete;

	static std::unique_ptr<microenvironment> create_from_config(const std::filesystem::path& config_file);

	void run_single_timestep() override;
	void serialize_state(real_t current_time) override;

	real_t get_substrate_density(index_t s, index_t x, index_t y, index_t z) const;

	void print_info(std::ostream& os) const;

	container_ptr agents;
	solver_ptr solver;
	serializer_ptr serializer;

	// environment configuration parameters
	std::string name, time_units, space_units;
	real_t diffusion_timestep;
	real_t simulation_time = 0.0;
	cartesian_mesh mesh;

	// diffusion-decay configuration parameters
	index_t substrates_count;
	std::vector<std::string> substrates_names;
	std::vector<std::string> substrates_units;
	std::unique_ptr<real_t[]> initial_conditions;
	std::unique_ptr<real_t[]> diffusion_coefficients;
	std::unique_ptr<real_t[]> decay_rates;

	// dirichlet interior configuration parameters
	index_t dirichlet_interior_voxels_count = 0;
	std::unique_ptr<index_t[]> dirichlet_interior_voxels;
	std::unique_ptr<real_t[]> dirichlet_interior_values;
	std::unique_ptr<bool[]> dirichlet_interior_conditions;

	// dirichlet boundary configuration parameters
	std::array<std::unique_ptr<real_t[]>, 3> dirichlet_min_boundary_values = { nullptr, nullptr, nullptr };
	std::array<std::unique_ptr<real_t[]>, 3> dirichlet_max_boundary_values = { nullptr, nullptr, nullptr };
	std::array<std::unique_ptr<bool[]>, 3> dirichlet_min_boundary_conditions = { nullptr, nullptr, nullptr };
	std::array<std::unique_ptr<bool[]>, 3> dirichlet_max_boundary_conditions = { nullptr, nullptr, nullptr };

	// bulk saturation-uptake configuration parameters
	std::unique_ptr<bulk_functor> bulk_fnc;

	// cell saturation-uptake configuration parameters
	bool compute_internalized_substrates;
};

} // namespace physicore::biofvm
