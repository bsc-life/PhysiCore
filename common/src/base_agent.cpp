#include "base_agent.h"

using namespace physicore;

base_agent::base_agent(index_t id, base_agent_data& data) : index(id), data(data) {}

std::span<real_t> base_agent::get_position()
{
	return std::span<real_t>(&data.positions[index * data.dims], data.dims);
}
