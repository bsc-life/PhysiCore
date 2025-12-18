#pragma once

#include <span>

#include "concepts.h"
#include "types.h"

namespace physicore {

template <derived_from_base_agent T>
class generic_agent_interface_container;

class base_agent_interface
{
protected:
	friend generic_agent_interface_container<base_agent_interface>;

	index_t index;

public:
	explicit base_agent_interface(index_t index) : index(index) {}

	virtual std::span<real_t> position() = 0;
	virtual ~base_agent_interface() = default;
};

} // namespace physicore
