#pragma once

#include <types.h>

#include <noarr/structures_extended.hpp>
#include <thrust/device_vector.h>

#include "../../../include/agent_data_generic_storage.h"
#include "../../../include/agent_generic_storage.h"
#include "../../../include/microenvironment.h"
#include "base_agent_data_generic_storage.h"
#include "base_agent_generic_storage.h"

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
	#define PHYSICORE_THRUST_DEVICE_FN __device__
	#define PHYSICORE_THRUST_STD cuda::std
#else
	#define PHYSICORE_THRUST_DEVICE_FN
	#define PHYSICORE_THRUST_STD std
#endif

namespace physicore::biofvm::kernels::thrust_solver {

using device_agent_data = agent_data_generic_storage<std::vector>;
using device_base_agent_data = base_agent_data_generic_storage<std::vector>;

using device_agent = agent_generic_storage<device_base_agent_data, device_agent_data>;
using device_base_agent = base_agent_generic_storage<device_base_agent_data>;

enum class data_residency
{
	HOST,
	DEVICE
};

class device_manager
{
	data_residency residency;

public:
	device_agent_data d_agent_data;
	std::unique_ptr<real_t> substrate_densities;

	void initialize(thrust::device_ptr<real_t> d_substrate_densities);

	void transfer_to_host(microenvironment& m);
	void transfer_to_device(microenvironment& m);
};


} // namespace physicore::biofvm::kernels::thrust_solver
