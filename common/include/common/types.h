#pragma once

#include <cstdint>

namespace physicore {

using real_t = double;
using index_t = std::uint64_t;
using sindex_t = std::int64_t;


#if defined(_MSC_VER)
	#define PHYSICORE_RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
	#define PHYSICORE_RESTRICT __restrict__
#else
	#define PHYSICORE_RESTRICT
#endif

} // namespace physicore
