#pragma once

#include <types.h>

#include <noarr/structures_extended.hpp>
#include <thrust/universal_vector.h>

#include "../../../include/agent_data.h"
#include "../../../include/microenvironment.h"

namespace physicore::biofvm::kernels::thrust_solver {

using device_agent_data = agent_data_generic<thrust::universal_vector>;

template <typename real_t>
using device_vector = thrust::universal_vector<real_t>;

struct device_data
{
	device_vector<real_t> substrate_densities;

	device_agent_data agent_data;

	auto get_substrates_layout(const microenvironment& m)
	{
		return noarr::scalar<real_t>()
			   ^ noarr::vectors<'s', 'x', 'y', 'z'>(m.substrates_count, m.mesh.grid_shape[0], m.mesh.grid_shape[1],
													m.mesh.grid_shape[2]);
	}
};

} // namespace physicore::biofvm::kernels::thrust_solver
