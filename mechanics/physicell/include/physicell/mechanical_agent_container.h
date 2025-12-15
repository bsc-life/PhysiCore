#pragma once

#include <common/base_agent.h>
#include <common/generic_agent_container.h>

#include "agent_data.h"
#include "mechanical_agent.h"

namespace physicore::mechanics::physicell {

/**
 * @brief Container for managing mechanical agents with synchronized data storage.
 *
 * Uses structure-of-arrays (SoA) pattern for efficient cache-friendly access to agent
 * properties (position, velocity, radius, adhesion, motility parameters, etc.).
 *
 * Agents are accessed as proxy objects (`base_agent` wrappers around `mechanical_agent_data`)
 * with automatic ID assignment on insertion.
 *
 * Thread-safe for reads; external synchronization needed for writes (add/remove).
 *
 * @see physicore::generic_agent_and_data_container for implementation details
 * @see physicore::biofvm::agent_container for similar pattern in reactions-diffusion module
 */
using mechanical_agent_container = physicore::generic_agent_and_data_container<physicore::base_agent, mechanical_agent>;

/**
 * @brief Interface-only container for mechanical agents (polymorphic access).
 *
 * Provides type-erased access to agents via `agent_interface`.
 * Useful for algorithms that need to work with agent interfaces without direct data access.
 */
using mechanical_agent_container_interface = physicore::generic_agent_interface_container<agent_interface>;

/**
 * @brief Shared pointer type for mechanical agent containers.
 */
using mechanical_container_ptr = std::shared_ptr<mechanical_agent_container_interface>;

} // namespace physicore::mechanics::physicell
