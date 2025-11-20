#pragma once

#include <common/base_agent_generic_storage.h>
#include <common/types.h>

#include "agent_interface.h"

namespace physicore::biofvm {

template <typename BaseAgentDataType, typename AgentDataType>
class agent_generic_storage : public physicore::base_agent_generic_storage<BaseAgentDataType>,
							  public virtual agent_interface
{
protected:
	AgentDataType& data;

public:
	using DataType = AgentDataType;
	using InterfaceType = agent_interface;

	agent_generic_storage(index_t index, AgentDataType& data)
		: base_agent_interface(index),
		  physicore::base_agent_generic_storage<BaseAgentDataType>(index, data.base_data),
		  data(data)
	{}

	agent_generic_storage(index_t index,
						  std::tuple<std::unique_ptr<BaseAgentDataType>, std::unique_ptr<AgentDataType>>& datas)
		: agent_generic_storage(index, *std::get<std::unique_ptr<AgentDataType>>(datas))
	{}

	std::span<real_t> secretion_rates() override
	{
		return std::span<real_t>(&data.secretion_rates[index * data.substrate_count], data.substrate_count);
	}

	std::span<real_t> saturation_densities() override
	{
		return std::span<real_t>(&data.saturation_densities[index * data.substrate_count], data.substrate_count);
	}

	std::span<real_t> uptake_rates() override
	{
		return std::span<real_t>(&data.uptake_rates[index * data.substrate_count], data.substrate_count);
	}

	std::span<real_t> net_export_rates() override
	{
		return std::span<real_t>(&data.net_export_rates[index * data.substrate_count], data.substrate_count);
	}

	std::span<real_t> internalized_substrates() override
	{
		return std::span<real_t>(&data.internalized_substrates[index * data.substrate_count], data.substrate_count);
	}

	std::span<real_t> fraction_released_at_death() override
	{
		return std::span<real_t>(&data.fraction_released_at_death[index * data.substrate_count], data.substrate_count);
	}

	std::span<real_t> fraction_transferred_when_ingested() override
	{
		return std::span<real_t>(&data.fraction_transferred_when_ingested[index * data.substrate_count],
								 data.substrate_count);
	}

	real_t& volume() override { return data.volumes[index]; }
};

} // namespace physicore::biofvm
