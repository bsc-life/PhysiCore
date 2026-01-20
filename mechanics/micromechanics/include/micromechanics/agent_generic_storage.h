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
	// Kinematics
	virtual std::span<real_t> velocity() = 0;
	virtual std::span<real_t> previous_velocity() = 0;

	// Geometry
	virtual real_t& radius() = 0;
	virtual std::uint8_t& is_movable() = 0;

	// Per-agent interaction strengths (can vary within a cell type)
	virtual real_t& cell_cell_adhesion_strength() = 0;
	virtual real_t& cell_cell_repulsion_strength() = 0;
	virtual real_t& relative_maximum_adhesion_distance() = 0;

	// Topology
	virtual std::span<index_t> neighbors() = 0;
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

	// Geometry
	real_t& radius() override { return data.radii[this->index]; }
	std::uint8_t& is_movable() override { return data.is_movable[this->index]; }

	// Per-agent interaction strengths
	real_t& cell_cell_adhesion_strength() override { return data.cell_cell_adhesion_strength[this->index]; }
	real_t& cell_cell_repulsion_strength() override { return data.cell_cell_repulsion_strength[this->index]; }
	real_t& relative_maximum_adhesion_distance() override
	{
		return data.relative_maximum_adhesion_distance[this->index];
	}

	// Topology
	std::span<index_t> neighbors() override { return std::span<index_t>(data.neighbors[this->index]); }
};

#ifdef _MSC_VER
	#pragma warning(pop)
#endif

} // namespace physicore::mechanics::micromechanics
