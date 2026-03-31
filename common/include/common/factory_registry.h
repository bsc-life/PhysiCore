#pragma once

#include <cassert>
#include <concepts>
#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace physicore {

/**
 * @brief Named factory registry for runtime-selectable backends.
 *
 * This is a small, generic utility used to implement BioFVM-style pluggable backends
 * (e.g., solvers). It is intentionally header-only and templated to avoid cross-module
 * link dependencies.
 */
template <typename ProductT>
class factory_registry
{
public:
	using product_ptr = std::unique_ptr<ProductT>;
	using factory_func_t = std::function<product_ptr()>;
	using registry_map_t = std::unordered_map<std::string, factory_func_t>;

	/**
	 * @brief Register a factory.
	 * @return true if name was newly registered, false if it already existed.
	 */
	bool register_factory(std::string name, factory_func_t&& factory)
	{
		auto [it, emplaced] = factories_.try_emplace(std::move(name), std::move(factory));
		return emplaced;
	}

	/**
	 * @brief Create an instance by name.
	 *
	 * In debug builds this asserts on unknown names (matching existing subsystem registries).
	 * In release builds it returns nullptr.
	 */
	product_ptr get(std::string_view name)
	{
		const auto it = factories_.find(std::string(name));
		if (it == factories_.end())
		{
			assert(false && "Unknown factory requested");
			return nullptr;
		}
		return it->second();
	}

	bool is_available(std::string_view name) const { return factories_.contains(std::string(name)); }

	std::vector<std::string> available_names() const
	{
		std::vector<std::string> names;
		names.reserve(factories_.size());
		for (const auto& [name, _] : factories_)
		{
			names.push_back(name);
		}
		return names;
	}

	// Compatibility with existing solver registries.
	std::vector<std::string> available_solvers() const { return available_names(); }

	static factory_registry& instance()
	{
		static factory_registry r;
		return r;
	}

private:
	registry_map_t factories_;
};

template <typename BaseT, typename ImplT>
concept derived_from_base = std::derived_from<ImplT, BaseT>;

/**
 * @brief Helper for static self-registration.
 */
template <typename BaseT, typename ImplT>
	requires derived_from_base<BaseT, ImplT>
struct factory_registry_adder
{
	explicit factory_registry_adder(std::string name)
	{
		factory_registry<BaseT>::instance().register_factory(std::move(name),
															 []() { return std::make_unique<ImplT>(); });
	}
};

} // namespace physicore
