#pragma once

#include <cassert>
#include <concepts>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

namespace physicore::common {
template <class Base>
class factory_registry
{
public:
	using ptr_t = std::unique_ptr<Base>;
	using factory_func_t = std::function<ptr_t()>;
	using map_t = std::unordered_map<std::string, factory_func_t>;

	bool register_factory(std::string name, factory_func_t&& f);
	ptr_t get(const std::string& name);

private:
	map_t factories_;
};

template <class Base>
bool factory_registry<Base>::register_factory(std::string name, factory_func_t&& f)
{
	auto [it, emplaced] = factories_.try_emplace(std::move(name), std::move(f));
	return emplaced;
}

template <class Base>
typename factory_registry<Base>::ptr_t factory_registry<Base>::get(const std::string& name)
{
	auto it = factories_.find(name);
	if (it == factories_.end())
	{
		assert(false);
		return nullptr;
	}

	return it->second();
}

template <class Derived, class Registry, class Base>
	requires std::derived_from<Derived, Base>
struct registry_adder
{
	explicit registry_adder(std::string name)
	{
		Registry::instance().register_factory(std::move(name), []() { return std::make_unique<Derived>(); });
	}
};

} // namespace physicore::common
