#pragma once

#include <cassert>
#include <concepts>
#include <cstring>
#include <vector>

#include "types.h"

namespace physicore {

template <typename T, typename ValueType>
concept VectorLike = requires(T t, ValueType val, std::size_t n) {
	typename T::value_type;
	{ t.push_back(val) };
	{ t.size() } -> std::convertible_to<std::size_t>;
	{ t[n] } -> std::convertible_to<ValueType&>;
	{ t.resize(n) };
};

template <template <typename...> typename ContainerType = std::vector>
struct base_agent_data_generic
{
	static_assert(VectorLike<ContainerType<real_t>, real_t>,
				  "ContainerType must satisfy VectorLike concept with real_t");

	index_t agents_count = 0;
	index_t dims;

	explicit base_agent_data_generic(index_t dims = 3) : dims(dims) {}

	void add()
	{
		++agents_count;
		positions.resize(agents_count * dims);
	}

	void remove_at(index_t position)
	{
		assert(position < agents_count);

		if (position >= agents_count)
			return;
		--agents_count;

		if (position != agents_count)
		{
			move_vector(&positions[position * dims], &positions[agents_count * dims], dims);
		}

		positions.resize(agents_count * dims);
	}

	template <typename T>
	static void move_scalar(T* dst, const T* src)
	{
		dst[0] = src[0];
	}

	template <typename T>
	static void move_vector(T* dst, const T* src, index_t size)
	{
		std::memcpy(dst, src, size * sizeof(T));
	}

	ContainerType<real_t> positions;
};

using base_agent_data = base_agent_data_generic<std::vector>;

} // namespace physicore
