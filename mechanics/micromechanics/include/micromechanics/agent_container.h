#pragma once

#include <memory>

#include <common/base_agent.h>
#include <common/generic_agent_container.h>

#include "agent.h"
#include "cell_interface.h"

namespace physicore::mechanics::micromechanics {

class agent_container : public generic_agent_and_data_container<base_agent, agent>
{
public:
	using generic_agent_and_data_container<base_agent, agent>::generic_agent_and_data_container;
};
using agent_container_interface = generic_agent_interface_container<cell_interface>;
using container_ptr = std::shared_ptr<agent_container_interface>;

} // namespace physicore::mechanics::micromechanics
