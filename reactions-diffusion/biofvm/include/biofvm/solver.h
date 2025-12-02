#pragma once

#include <memory>

#include <biofvm/biofvm_export.h>
#include <common/types.h>

namespace physicore::biofvm {

class microenvironment;

class BIOFVM_EXPORT solver
{
public:
	// Set initial values (such as substrate densities) in the microenvironment
	virtual void initialize(microenvironment& m) = 0;

	// Solve the diffusion-decay equations for a given number of iterations
	virtual void solve(microenvironment& m, index_t iterations) = 0;

	// Get the substrate density at a given voxel
	virtual real_t get_substrate_density(index_t s, index_t x, index_t y, index_t z) const = 0;
	virtual real_t& get_substrate_density(index_t s, index_t x, index_t y, index_t z) = 0;

	// Transfer data to/from device (if applicable)
	virtual void transfer_to_device([[maybe_unused]] microenvironment& m) { /* Default host solver */ }
	virtual void transfer_to_host([[maybe_unused]] microenvironment& m) { /* Default host solver */ }

	// Reinitialize Dirichlet conditions (e.g., after modifying boundary or interior conditions)
	virtual void reinitialize_dirichlet(microenvironment& m) = 0;

	virtual ~solver() = default;
};

using solver_ptr = std::unique_ptr<solver>;

} // namespace physicore::biofvm
