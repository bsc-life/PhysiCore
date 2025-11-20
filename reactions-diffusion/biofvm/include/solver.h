#pragma once

#include <memory>

#include "types.h"

namespace physicore::biofvm {

class microenvironment;

class solver
{
public:
	// Set initial values (such as substrate densities) in the microenvironment
	virtual void initialize(microenvironment& m) = 0;

	// Solve the diffusion-decay equations for a given number of iterations
	virtual void solve(microenvironment& m, index_t iterations) = 0;

	// Get the substrate density at a given voxel
	virtual real_t get_substrate_density(index_t s, index_t x, index_t y, index_t z) const = 0;

	virtual ~solver() = default;
};

using solver_ptr = std::unique_ptr<solver>;

} // namespace physicore::biofvm
