#pragma once

#include <common/base_agent_data.h>

#include "agent_data.h"
#include "cell_data.h"

namespace physicore::mechanics::micromechanics {

/// Reset per-cell aggregate buffers that are computed from agents.
void reset_cell_aggregates(cell_data& cells);

/// Aggregate COM position/velocity, per-compartment counts, and a pressure proxy from agent forces.
void aggregate_cell_data_from_agents(const physicore::base_agent_data& base,
								 const agent_data& agents,
								 cell_data& cells);

} // namespace physicore::mechanics::micromechanics

