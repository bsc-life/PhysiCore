#include "environment.h"

#include <memory>
#include <stdexcept>
#include <tuple>

#include <common/base_agent_data.h>

using namespace physicore::mechanics::physicell;

environment::environment(real_t timestep, index_t dims, index_t agent_types_count, index_t substrates_count)
	: timestep(timestep)
{
	auto base_data = std::make_unique<physicore::base_agent_data>(dims);
	auto data = std::make_unique<mechanical_agent_data>(*base_data, agent_types_count, substrates_count);
	agents = std::make_unique<mechanical_agent_container>(std::move(base_data), std::move(data));
}

void environment::run_single_timestep()
{
	if (solver)
	{
		solver->solve(*this, 1);
		return;
	}

	(void)timestep;
}

void environment::serialize_state(real_t current_time)
{
	if (serializer)
	{
		serializer->serialize(current_time);
	}
}

mechanical_agent_data& environment::get_agent_data()
{
	if (!agents)
	{
		throw std::runtime_error("environment has no mechanical agent container");
	}
	auto& data_ptr = std::get<std::unique_ptr<mechanical_agent_data>>(agents->agent_datas);
	return *data_ptr;
}
