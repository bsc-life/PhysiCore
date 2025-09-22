#pragma once

#include <vector>

#include "../../../../include/microenvironment.h"

namespace physicore::biofvm::kernels::cpu {

struct problem_t
{
	index_t dims = 3;
	index_t dx = 20, dy = 20, dz = 20;
	index_t nx = 1, ny = 1, nz = 1;
	index_t substrates_count = 1;

	index_t iterations = 1;

	real_t dt = 0.01;
	std::vector<real_t> diffusion_coefficients;
	std::vector<real_t> decay_rates;
	std::vector<real_t> initial_conditions;

	static problem_t construct(const microenvironment& m, index_t iterations)
	{
		problem_t problem;
		problem.dims = m.mesh.dims;
		problem.dx = m.mesh.voxel_shape[0];
		problem.dy = m.mesh.voxel_shape[1];
		problem.dz = m.mesh.voxel_shape[2];
		problem.nx = m.mesh.grid_shape[0];
		problem.ny = m.mesh.grid_shape[1];
		problem.nz = m.mesh.grid_shape[2];
		problem.substrates_count = m.substrates_count;
		problem.iterations = iterations;
		problem.dt = static_cast<real_t>(m.diffusion_timestep);
		problem.diffusion_coefficients =
			std::vector<real_t>(m.diffusion_coefficients.get(), m.diffusion_coefficients.get() + m.substrates_count);
		problem.decay_rates = std::vector<real_t>(m.decay_rates.get(), m.decay_rates.get() + m.substrates_count);
		problem.initial_conditions =
			std::vector<real_t>(m.initial_conditions.get(), m.initial_conditions.get() + m.substrates_count);
		return problem;
	}
};

} // namespace physicore::biofvm::kernels::cpu
