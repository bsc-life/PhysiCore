#include "base_agent_data.h"

#include <cassert>

using namespace physicore;

base_agent_data::base_agent_data(index_t dims) : dims(dims) {}

void base_agent_data::add()
{
	++agents_count;
	positions.resize(agents_count * dims);
}

void base_agent_data::remove_at(index_t position)
{
	assert(position < agents_count);

	if (position >= agents_count)
		return;
	--agents_count;

	if (position != agents_count)
	{
		move_vector(&positions[position * dims], &positions[agents_count * dims], dims);
	}

	positions.resize(agents_count * dims);
}
