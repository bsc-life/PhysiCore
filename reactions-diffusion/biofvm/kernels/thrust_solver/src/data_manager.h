#pragma once

#include <biofvm/agent.h>
#include <biofvm/agent_container.h>
#include <biofvm/microenvironment.h>
#include <common/generic_agent_container.h>
#include <common/generic_agent_solver.h>
#include <common/types.h>
#include <noarr/structures_extended.hpp>
#include <thrust/device_vector.h>

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
	#define PHYSICORE_THRUST_DEVICE_FN __device__
	#define PHYSICORE_THRUST_STD cuda::std
#else
	#define PHYSICORE_THRUST_DEVICE_FN
	#define PHYSICORE_THRUST_STD std
#endif

namespace physicore::biofvm::kernels::thrust_solver {

class diffusion_solver;

enum class data_residency
{
	HOST,
	DEVICE
};

// This class takes care of moving data between host and device memory.
// As of now, substrate_densities default residency is device while agent_data default residency is host.
// Therefore, the exposed public pointers to their residency conterparts. I.e., host densities and device agent_data.
class data_manager : private generic_agent_solver<agent>
{
	data_residency residency;

	std::size_t densities_size_bytes;

	thrust::device_vector<real_t> d_positions;
	thrust::device_vector<real_t> d_secretion_rates;
	thrust::device_vector<real_t> d_saturation_densities;
	thrust::device_vector<real_t> d_uptake_rates;
	thrust::device_vector<real_t> d_net_export_rates;
	thrust::device_vector<real_t> d_internalized_substrates;
	thrust::device_vector<real_t> d_fraction_released_at_death;
	thrust::device_vector<real_t> d_fraction_transferred_when_ingested;
	thrust::device_vector<real_t> d_volumes;
	thrust::device_ptr<real_t> d_substrate_densities;

	std::shared_ptr<agent_container_interface> h_agent_container;
	std::unique_ptr<real_t[]> h_substrate_densities;

public:
	void initialize(const microenvironment& m, diffusion_solver& d_solver);

	void transfer_to_host();
	void transfer_to_device();

	real_t* substrate_densities = nullptr;

	real_t* positions = nullptr;
	real_t* secretion_rates = nullptr;
	real_t* saturation_densities = nullptr;
	real_t* uptake_rates = nullptr;
	real_t* net_export_rates = nullptr;
	real_t* internalized_substrates = nullptr;
	real_t* fraction_released_at_death = nullptr;
	real_t* fraction_transferred_when_ingested = nullptr;
	real_t* volumes = nullptr;
};


} // namespace physicore::biofvm::kernels::thrust_solver
