#pragma once

#include "data_manager.h"

namespace physicore::biofvm::kernels::thrust_solver {

struct device_bulk_functor
{
	PHYSICORE_THRUST_DEVICE_FN virtual real_t supply_rates(index_t s, index_t x, index_t y, index_t z) = 0;
	PHYSICORE_THRUST_DEVICE_FN virtual real_t uptake_rates(index_t s, index_t x, index_t y, index_t z) = 0;
	PHYSICORE_THRUST_DEVICE_FN virtual real_t supply_target_densities(index_t s, index_t x, index_t y, index_t z) = 0;
	PHYSICORE_THRUST_DEVICE_FN virtual ~device_bulk_functor() = default;
};

} // namespace physicore::biofvm::kernels::thrust_solver
