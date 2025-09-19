#include "base_agent.h"

#include <cassert>

using namespace physicore;

base_agent::base_agent(index_t id, base_agent_data& data) : index(id), data(data) { assert(index < data.agents_count); }

std::span<real_t> base_agent::get_position()
{
	return std::span<real_t>(&data.positions[index * data.dims], data.dims);
}
