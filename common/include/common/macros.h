#pragma once

namespace physicore {

#if defined(_MSC_VER)
	#define PHYSICORE_RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
	#define PHYSICORE_RESTRICT __restrict__
#else
	#define PHYSICORE_RESTRICT
#endif

} // namespace physicore
