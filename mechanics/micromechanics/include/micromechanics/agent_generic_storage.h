#pragma once

#include <cstdint>
#include <span>

#include <common/base_agent_generic_storage.h>
#include <common/types.h>

#include "agent_data.h"

namespace physicore::mechanics::micromechanics {

namespace internal {

/**
 * @brief Internal interface for mechanics agents (sub-cellular compartments).
 *
 * This is intentionally not part of the public micromechanics API surface.
 */
class agent_interface : public virtual base_agent_interface
{
public:
	// Agent Classification
	virtual std::uint8_t& compartment_type() = 0;
	virtual index_t& cell_id() = 0;

	// Kinematics
	virtual std::span<real_t> velocity() = 0;
	virtual std::span<real_t> previous_velocity() = 0;
	virtual std::span<real_t> force() = 0;

	// Topology (Kelvin-Voigt)
	virtual std::span<index_t> spring_attachments() = 0;
};

} // namespace internal

#ifdef _MSC_VER
	#pragma warning(push)
	#pragma warning(disable : 4250) // inherits via dominance (intentional design)
#endif

/**
 * @brief Template implementation of agent storage for mechanics agents.
 *
 * Provides access to agent-level mechanics properties stored in SoA format.
 * Cell-level properties (pressure, volume, etc.) are accessed via cell_data.
 */
template <typename BaseAgentDataType, typename AgentDataType>
class agent_generic_storage : public physicore::base_agent_generic_storage<BaseAgentDataType>,
							  public virtual internal::agent_interface
{
protected:
	AgentDataType& data;

public:
	using DataType = AgentDataType;
	using InterfaceType = internal::agent_interface;

	agent_generic_storage(index_t index, AgentDataType& data)
		: base_agent_interface(index),
		  physicore::base_agent_generic_storage<BaseAgentDataType>(index, data.base_data),
		  data(data)
	{}

	agent_generic_storage(index_t index,
						  std::tuple<std::unique_ptr<BaseAgentDataType>, std::unique_ptr<AgentDataType>>& datas)
		: agent_generic_storage(index, *std::get<std::unique_ptr<AgentDataType>>(datas))
	{}

	// Kinematics
	std::uint8_t& compartment_type() override { return data.compartment_types[this->index]; }
	index_t& cell_id() override { return data.cell_ids[this->index]; }

	std::span<real_t> velocity() override
	{
		const index_t dims = data.base_data.dims;
		return std::span<real_t>(&data.velocities[this->index * dims], dims);
	}

	std::span<real_t> previous_velocity() override
	{
		const index_t dims = data.base_data.dims;
		return std::span<real_t>(&data.previous_velocities[this->index * dims], dims);
	}

	std::span<real_t> force() override
	{
		const index_t dims = data.base_data.dims;
		return std::span<real_t>(&data.forces[this->index * dims], dims);
	}

	std::span<index_t> spring_attachments() override
	{
		return std::span<index_t>(data.spring_attachments[this->index]);
	}
};

#ifdef _MSC_VER
	#pragma warning(pop)
#endif

} // namespace physicore::mechanics::micromechanics
