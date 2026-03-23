#pragma once

#include "types.h"

namespace physicore {

class random
{
public:
	static random& instance();

	static real_t uniform(const real_t min = 0, const real_t max = 1);

	static real_t normal(const real_t mean = 0, const real_t std = 1);

	static void set_seed(unsigned int seed);
};

} // namespace physicore
