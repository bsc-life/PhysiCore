#pragma once

#include <vector>

#include "types.h"

namespace physicore {
struct base_agent_data
{
	index_t agents_count = 0;
	index_t dims = 3; // Default to 3D

	void add();
	void remove_at(index_t position);

	std::vector<real_t> positions;
};
} // namespace physicore
