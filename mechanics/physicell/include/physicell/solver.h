#pragma once

#include <memory>

#include <common/types.h>
#include <physicell/physicell_export.h>

namespace physicore::mechanics::physicell {

class environment;

class PHYSICELL_EXPORT solver
{
public:
	virtual void initialize(environment& e) = 0;

	virtual void solve(environment& e, index_t iterations) = 0;

	virtual ~solver() = default;
};

using solver_ptr = std::unique_ptr<solver>;

} // namespace physicore::mechanics::physicell
