#pragma once

#include <biofvm/microenvironment.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>

#include "diffusion_solver.h"
#include "namespace_config.h"

/*
Performs secretion and uptake of cells.
Updates substrate denisities of the cell's voxel and conditionally updates the cell's internalized substrates.

D = substrate densities
I = internalized substrates
S = secretion rates
U = uptake rates
T = saturation densities
N = net export rates
c = cell_volume
v = voxel_volume

I -= v((-(c/v)*dt*(U+S)*D + (c/v)*dt*S*T) / (1 + (c/v)*dt*(U+S)) + (1/v)dt*N)

Updating substrate densities is more complex, one has to take care about the case when we have more cells in the same
voxel. The following formula is used:

D = (D + sum_k{(c_k/v)*dt*S_k*T_k}) / (1 + sum_k{(c_k/v)*dt*(U_k+S_k)}) + sum_k{(1/v)*dt*N_k}

where sum is over all cells in the voxel.

Also handles release of internalized substrates:

F = fraction released at death

D = D + I*F/v
*/

namespace physicore::biofvm::kernels::PHYSICORE_THRUST_SOLVER_NAMESPACE {

class cell_solver
{
	bool compute_internalized_substrates_;

	thrust::device_vector<real_t> numerators_;
	thrust::device_vector<real_t> denominators_;
	thrust::device_vector<real_t> factors_;

	thrust::device_vector<real_t> reduced_numerators_;
	thrust::device_vector<real_t> reduced_denominators_;
	thrust::device_vector<real_t> reduced_factors_;

	thrust::device_ptr<bool> is_conflict_;

	thrust::device_ptr<index_t> ballots_;

	void resize(const microenvironment& m);

public:
	cell_solver() = default;
	cell_solver(const cell_solver&) = delete;
	cell_solver& operator=(const cell_solver&) = delete;

	void initialize(const microenvironment& m);

	void simulate_secretion_and_uptake(microenvironment& m, diffusion_solver& d_solver, data_manager& data,
									   bool recompute);

	void release_internalized_substrates(const microenvironment& m, diffusion_solver& d_solver, data_manager& data,
										 index_t agent_index);

	void release_internalized_substrates(const microenvironment& m, diffusion_solver& d_solver, data_manager& data);

	~cell_solver();
};

} // namespace physicore::biofvm::kernels::PHYSICORE_THRUST_SOLVER_NAMESPACE
