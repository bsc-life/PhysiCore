#pragma once

#include <concepts>

#include <common/types.h>

namespace physicore {

class base_agent_interface;

template <typename T>
concept agent_data_type = requires(T data, index_t pos) {
	{ data.add() } -> std::same_as<void>;
	{ data.remove_at(pos) } -> std::same_as<void>;
};

// Concept to check if a type is derived from base_agent
template <typename T>
concept derived_from_base_agent = std::derived_from<T, base_agent_interface>;

} // namespace physicore
