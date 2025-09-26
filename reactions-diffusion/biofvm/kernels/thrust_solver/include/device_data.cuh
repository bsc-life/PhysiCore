#pragma once

#include <types.h>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/universal_vector.h>

#include "../../../include/agent_data.h"

namespace physicore::biofvm::kernels::thrust_solver {

using device_agent_data = agent_data_generic<thrust::device_vector>;

struct data_manager
{
	std::unique_ptr<real_t[]> substrate_densities;
	thrust::universal_vector<real_t> d_substrate_densities;

	device_agent_data agent_data;
};

} // namespace physicore::biofvm::kernels::thrust_solver
