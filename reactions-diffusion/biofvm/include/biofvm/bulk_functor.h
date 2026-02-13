#pragma once

#include <biofvm/biofvm_export.h>
#include <common/types.h>

namespace physicore::reactions_diffusion::biofvm {

struct BIOFVM_EXPORT bulk_functor
{
	virtual real_t supply_rates(index_t s, index_t x, index_t y, index_t z) = 0;
	virtual real_t uptake_rates(index_t s, index_t x, index_t y, index_t z) = 0;
	virtual real_t supply_target_densities(index_t s, index_t x, index_t y, index_t z) = 0;
	virtual ~bulk_functor() = default;
};

} // namespace physicore::reactions_diffusion::biofvm
