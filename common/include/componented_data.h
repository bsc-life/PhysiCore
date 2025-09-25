#pragma once

#include <cassert>

#include "data_component.h"

namespace physicore {

template <typename T>
concept derived_from_data_component = std::derived_from<T, data_component>;

struct componented_data
{
	component_map_t components;

	template <derived_from_data_component... ComponentTypes>
	componented_data([[maybe_unused]] ComponentTypes...)
	{
		(add_component<ComponentTypes>(), ...);
	}

	template <derived_from_data_component ComponentType>
	ComponentType* get()
	{
		auto it = components.find(typeid(ComponentType));
		if (it == components.end())
		{
			assert(false);
			return nullptr;
		}

		return static_cast<ComponentType*>(it->second.get());
	}

	template <derived_from_data_component ComponentType>
	bool add_component()
	{
		auto [it, emplaced] = components.try_emplace(typeid(ComponentType), std::make_unique<ComponentType>());

		return emplaced;
	}
};

} // namespace physicore
