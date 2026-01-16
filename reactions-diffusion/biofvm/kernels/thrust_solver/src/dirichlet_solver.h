#pragma once

#include <thrust/device_vector.h>

#include "diffusion_solver.h"
#include "namespace_config.h"

/*
This solver applies Dirichlet boundary conditions to the microenvironment.
I.e., it sets the values of the substrate concentrations at specific voxels to a constant value.

Implementation works with 3 arrays:
m.dirichlet_voxels - array of dirichlet voxel (1D/2D/3D) indices
m.dirichlet_conditions - array of bools specifying if a substrate of a dirichlet voxel has a dirichled codition
m.dirichlet_values - array of dirichlet values for each substrate with a dirichlet condition
*/

namespace physicore::reactions_diffusion::biofvm::kernels::PHYSICORE_THRUST_SOLVER_NAMESPACE {

class dirichlet_solver
{
	thrust::device_vector<index_t> dirichlet_interior_voxels;
	thrust::device_vector<real_t> dirichlet_interior_values;
	thrust::device_vector<bool> dirichlet_interior_conditions;

	std::array<thrust::device_vector<real_t>, 3> dirichlet_min_boundary_values;
	std::array<thrust::device_vector<real_t>, 3> dirichlet_max_boundary_values;
	std::array<thrust::device_vector<bool>, 3> dirichlet_min_boundary_conditions;
	std::array<thrust::device_vector<bool>, 3> dirichlet_max_boundary_conditions;

public:
	void initialize(microenvironment& m);
	void solve(microenvironment& m, diffusion_solver& d_solver);
};

} // namespace physicore::reactions_diffusion::biofvm::kernels::PHYSICORE_THRUST_SOLVER_NAMESPACE
