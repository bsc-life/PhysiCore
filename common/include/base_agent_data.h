#pragma once

#include <vector>

#include "types.h"

namespace physicore {
struct base_agent_data
{
	index_t agents_count;
	index_t dims;

	void add();
	void remove_at(index_t position);

	std::vector<real_t> positions;
};
} // namespace physicore