#include "physicell/solver_registry.h"

namespace physicore::mechanics::physicell {

solver_registry& solver_registry::instance()
{
	static solver_registry r;
	return r;
}

} // namespace physicore::mechanics::physicell
