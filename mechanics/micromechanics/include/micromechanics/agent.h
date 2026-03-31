#pragma once
#include <common/base_agent_data.h>

#include "agent_data.h"
#include "agent_generic_storage.h"

namespace physicore::mechanics::micromechanics {

using agent = agent_generic_storage<base_agent_data, agent_data>;

} // namespace physicore::mechanics::micromechanics
