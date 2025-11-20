#pragma once

#include "diffusion_solver.h"

/*
This solver applies Dirichlet boundary conditions to the microenvironment.
I.e., it sets the values of the substrate concentrations at specific voxels to a constant value.

Implementation works with 3 arrays:
m.dirichlet_voxels - array of dirichlet voxel (1D/2D/3D) indices
m.dirichlet_conditions - array of bools specifying if a substrate of a dirichlet voxel has a dirichled codition
m.dirichlet_values - array of dirichlet values for each substrate with a dirichlet condition
*/

namespace physicore::biofvm::kernels::openmp_solver {

class dirichlet_solver
{
public:
	static void solve(microenvironment& m, diffusion_solver& d_solver);
};

} // namespace physicore::biofvm::kernels::openmp_solver
