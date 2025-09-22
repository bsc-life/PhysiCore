#pragma once

#include "types.h"

namespace physicore::biofvm {

class microenvironment;

class solver
{
public:
	virtual void solve(microenvironment& m, index_t iterations) = 0;

	virtual ~solver() = default;
};

} // namespace physicore::biofvm
