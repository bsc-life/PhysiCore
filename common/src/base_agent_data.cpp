#include "base_agent_data.h"

using namespace physicore;

void base_agent_data::add()
{
	++agents_count;
	positions.resize(agents_count * dims);
}

void base_agent_data::remove_at(index_t position)
{
	if (position >= agents_count)
		return;
	--agents_count;

	for (size_t d = 0; d < dims; ++d)
	{
		positions[position * dims + d] = positions[agents_count * dims + d];
	}

	positions.resize(agents_count * dims);
}
