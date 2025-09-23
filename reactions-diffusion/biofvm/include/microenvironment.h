#pragma once

#include <functional>
#include <process.h>
#include <string>
#include <types.h>
#include <vector>

#include "agent_container.h"
#include "mesh.h"

namespace physicore::biofvm {

class microenvironment : public timestep_executor
{
public:
	microenvironment(const cartesian_mesh& mesh, index_t substrates_count, real_t timestep);

	microenvironment(microenvironment&&) = delete;
	microenvironment(const microenvironment&) = delete;
	microenvironment& operator=(const microenvironment&) = delete;
	microenvironment& operator=(microenvironment&&) = delete;

	void run_single_timestep() override;

	std::string name, time_units, space_units;
	std::vector<std::string> substrates_names;
	std::vector<std::string> substrates_units;

	cartesian_mesh mesh;

	agent_container agents;

	index_t substrates_count;
	real_t diffusion_timestep;

	std::unique_ptr<real_t[]> initial_conditions;
	std::unique_ptr<real_t[]> diffusion_coefficients;
	std::unique_ptr<real_t[]> decay_rates;

	index_t dirichlet_interior_voxels_count = 0;
	std::unique_ptr<index_t[]> dirichlet_interior_voxels;
	std::unique_ptr<real_t[]> dirichlet_interior_values;
	std::unique_ptr<bool[]> dirichlet_interior_conditions;

	std::array<std::unique_ptr<real_t[]>, 3> dirichlet_min_boundary_values = { nullptr, nullptr, nullptr };
	std::array<std::unique_ptr<real_t[]>, 3> dirichlet_max_boundary_values = { nullptr, nullptr, nullptr };
	std::array<std::unique_ptr<bool[]>, 3> dirichlet_min_boundary_conditions = { nullptr, nullptr, nullptr };
	std::array<std::unique_ptr<bool[]>, 3> dirichlet_max_boundary_conditions = { nullptr, nullptr, nullptr };

	using bulk_func_t = std::function<void(microenvironment& m, std::array<index_t, 3> voxel_idx, real_t* out)>;
	bulk_func_t supply_rate_func, uptake_rate_func, supply_target_densities_func;

	bool compute_internalized_substrates;
};

} // namespace physicore::biofvm
