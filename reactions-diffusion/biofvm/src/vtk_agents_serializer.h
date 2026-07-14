#pragma once

#include <string>
#include <string_view>
#include <vector>

#include <common/generic_agent_solver.h>
#include <common/types.h>

#include "agent.h"
#include "serializer.h"
#include "vtk_serializer_base.h"

namespace physicore::biofvm {

class vtk_agents_serializer : public vtk_serializer_base, public serializer, private generic_agent_solver<agent>
{
	index_t substrate_count_;
	std::vector<std::string> substrate_names_;

public:
	vtk_agents_serializer(std::string_view output_dir, const microenvironment& m);

	void serialize(const microenvironment& m, real_t current_time) override;
};

} // namespace physicore::biofvm
