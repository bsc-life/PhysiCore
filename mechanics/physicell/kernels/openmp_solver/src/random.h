#pragma once

#include <common/types.h>

namespace physicore {

class random
{
public:
	static real_t uniform(real_t min = 0, real_t max = 1);

	static real_t normal(real_t mean = 0, real_t std = 1);

	static void set_seed(unsigned int seed);
};

} // namespace physicore
