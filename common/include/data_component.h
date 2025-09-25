#pragma once

#include <memory>
#include <typeindex>
#include <unordered_map>

#include "types.h"

namespace physicore {
struct data_component
{
	virtual void resize(index_t n) = 0;
	virtual void remove_at(index_t position) = 0;

	virtual ~data_component() = default;
};

using component_map_t = std::unordered_map<std::type_index, std::unique_ptr<data_component>>;
} // namespace physicore
