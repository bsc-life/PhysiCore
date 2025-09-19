#pragma once

#include <cstring>
#include <vector>

#include "types.h"

namespace physicore {
struct base_agent_data
{
	index_t agents_count = 0;
	index_t dims = 3; // Default to 3D

	void add();
	void remove_at(index_t position);

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

	std::vector<real_t> positions;
};
} // namespace physicore
