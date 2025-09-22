#include "base_agent.h"

#include <cassert>

using namespace physicore;

base_agent::base_agent(index_t id, base_agent_data& data) : index(id), base_data(data)
{
	assert(index < data.agents_count);
}

std::span<real_t> base_agent::position()
{
	return std::span<real_t>(&base_data.positions[index * base_data.dims], base_data.dims);
}
