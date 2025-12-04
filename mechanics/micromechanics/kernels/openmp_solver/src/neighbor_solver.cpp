#include "neighbor_solver.h"

#include <micromechanics/environment.h>
#include <micromechanics/spatial_index.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void neighbor_solver::initialize(environment& /*e*/)
{
	if (initialized_)
		return;

	initialized_ = true;
}

void neighbor_solver::update_neighbors(environment& e)
{
	// Build the spatial index for efficient neighbor queries
	if (e.index)
	{
		e.index->build(e);
	}

	// Note: The actual neighbor lists are built on-demand during force calculation
	// using the spatial index. This is more efficient than pre-computing all neighbors
	// because different potentials may have different interaction distances.
	//
	// For Kelvin-Voigt model, we could pre-compute topology here if needed.
	// The legacy code had a static topology flag (cells_updated) but we use
	// dynamic neighbor finding for flexibility.
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
