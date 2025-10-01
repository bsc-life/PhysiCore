#pragma once

#include <types.h>

#include <noarr/structures_extended.hpp>
#include <thrust/universal_vector.h>

#include "../../../include/agent_data_generic_storage.h"
#include "../../../include/agent_generic_storage.h"
#include "base_agent_data_generic_storage.h"
#include "base_agent_generic_storage.h"

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
	#define PHYSICORE_DEVICE __device__
#else
	#define PHYSICORE_DEVICE
#endif

namespace physicore::biofvm::kernels::thrust_solver {

using device_agent_data = agent_data_generic_storage<thrust::universal_vector>;
using device_base_agent_data = base_agent_data_generic_storage<thrust::universal_vector>;

using device_agent = agent_generic_storage<device_base_agent_data, device_agent_data>;
using device_base_agent = base_agent_generic_storage<device_base_agent_data>;

template <typename real_t>
using device_vector = thrust::universal_vector<real_t>;


} // namespace physicore::biofvm::kernels::thrust_solver
