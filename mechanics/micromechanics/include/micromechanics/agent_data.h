#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

#include <common/base_agent_data.h>
#include <common/types.h>

namespace physicore::mechanics::micromechanics {

template <template <typename...> typename ContainerType = std::vector>
struct agent_data_generic_storage
{
public:
	physicore::base_agent_data_generic_storage<ContainerType>& base_data;

	// Agent Classification
	ContainerType<std::uint8_t> compartment_types; // Agent type for type-based interactions (0-255)
	ContainerType<index_t> cell_ids;			   // Cell ID this agent belongs to (-1 if standalone)

	// Physics State
	ContainerType<real_t> velocities;		   // n * dims
	ContainerType<real_t> previous_velocities; // n * dims (for Adams-Bashforth)
	ContainerType<real_t> forces; // n * dims magnitude of forces acting on agent (we can keep it for later)


	// Topology (Kelvin-Voigt)
	// Note: Nested vectors are not ideal for SoA but required for dynamic topology

	// Spring attachments (list of attached agent indices per agent)
	ContainerType<std::vector<index_t>> spring_attachments;

	index_t agents_count = 0;

	explicit agent_data_generic_storage(physicore::base_agent_data_generic_storage<ContainerType>& base_data);

	void add();
	void remove_at(index_t position);
};

template <template <typename...> typename ContainerType>
agent_data_generic_storage<ContainerType>::agent_data_generic_storage(
	physicore::base_agent_data_generic_storage<ContainerType>& base_data)
	: base_data(base_data)
{}

template <template <typename...> typename ContainerType>
void agent_data_generic_storage<ContainerType>::add()
{
	++agents_count;
	index_t dims = base_data.dims;

	// Agent Classification
	compartment_types.resize(agents_count, 0);
	cell_ids.resize(agents_count, static_cast<index_t>(-1)); // -1 = standalone agent

	velocities.resize(agents_count * dims, 0.0);
	previous_velocities.resize(agents_count * dims, 0.0);
	forces.resize(agents_count * dims, 0.0);

	spring_attachments.resize(agents_count);
}

template <template <typename...> typename ContainerType>
void agent_data_generic_storage<ContainerType>::remove_at(index_t position)
{
	assert(position < agents_count);
	if (position >= agents_count)
		return;

	--agents_count;
	index_t dims = base_data.dims;
	using base_move = physicore::base_agent_data_generic_storage<ContainerType>;

	if (position < agents_count)
	{
		// Agent Classification
		base_move::move_scalar(&compartment_types[position], &compartment_types[agents_count]);
		base_move::move_scalar(&cell_ids[position], &cell_ids[agents_count]);

		base_move::move_vector(&velocities[position * dims], &velocities[agents_count * dims], dims);
		base_move::move_vector(&previous_velocities[position * dims], &previous_velocities[agents_count * dims], dims);
		base_move::move_vector(&forces[position * dims], &forces[agents_count * dims], dims);

		spring_attachments[position] = std::move(spring_attachments[agents_count]);
	}

	// Resize to shrink
	compartment_types.resize(agents_count);
	cell_ids.resize(agents_count);
	velocities.resize(agents_count * dims);
	previous_velocities.resize(agents_count * dims);
	forces.resize(agents_count * dims);

	spring_attachments.resize(agents_count);
}

using agent_data = agent_data_generic_storage<std::vector>;

} // namespace physicore::mechanics::micromechanics
