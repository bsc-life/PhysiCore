#include "data_manager.h"

#include <memory>

#include "diffusion_solver.h"
#include "namespace_config.h"

using namespace physicore;
using namespace physicore::biofvm::kernels::PHYSICORE_THRUST_SOLVER_NAMESPACE;

void data_manager::initialize(const microenvironment& m, diffusion_solver& d_solver)
{
	h_agent_container = m.agents;

	densities_size_bytes = d_solver.get_substrates_layout() | noarr::get_size();
	d_substrate_densities = d_solver.get_substrates_pointer();

	// move densities to host such that during the initialization, all data is in one memory space - host
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
	{
		h_substrate_densities = std::make_unique<real_t[]>(densities_size_bytes / sizeof(real_t));

		thrust::copy_n(d_substrate_densities, densities_size_bytes / sizeof(real_t), h_substrate_densities.get());

		substrate_densities = h_substrate_densities.get();
	}
#else
	substrate_densities = d_substrate_densities.get();
#endif

	residency = data_residency::HOST;
}

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
void data_manager::transfer_to_host()
{
	auto& h_agent_data = retrieve_agent_data(*h_agent_container);

	thrust::copy(d_positions.begin(), d_positions.end(), h_agent_data.base_data.positions.begin());
	thrust::copy(d_secretion_rates.begin(), d_secretion_rates.end(), h_agent_data.secretion_rates.begin());
	thrust::copy(d_saturation_densities.begin(), d_saturation_densities.end(),
				 h_agent_data.saturation_densities.begin());
	thrust::copy(d_uptake_rates.begin(), d_uptake_rates.end(), h_agent_data.uptake_rates.begin());
	thrust::copy(d_net_export_rates.begin(), d_net_export_rates.end(), h_agent_data.net_export_rates.begin());
	thrust::copy(d_internalized_substrates.begin(), d_internalized_substrates.end(),
				 h_agent_data.internalized_substrates.begin());
	thrust::copy(d_fraction_released_at_death.begin(), d_fraction_released_at_death.end(),
				 h_agent_data.fraction_released_at_death.begin());
	thrust::copy(d_fraction_transferred_when_ingested.begin(), d_fraction_transferred_when_ingested.end(),
				 h_agent_data.fraction_transferred_when_ingested.begin());
	thrust::copy(d_volumes.begin(), d_volumes.end(), h_agent_data.volumes.begin());

	thrust::copy_n(d_substrate_densities, densities_size_bytes / sizeof(real_t), h_substrate_densities.get());

	substrate_densities = h_substrate_densities.get();

	residency = data_residency::HOST;
}

void data_manager::transfer_to_device()
{
	auto& h_agent_data = retrieve_agent_data(*h_agent_container);

	d_positions = h_agent_data.base_data.positions;
	d_secretion_rates = h_agent_data.secretion_rates;
	d_saturation_densities = h_agent_data.saturation_densities;
	d_uptake_rates = h_agent_data.uptake_rates;
	d_net_export_rates = h_agent_data.net_export_rates;
	d_internalized_substrates = h_agent_data.internalized_substrates;
	d_fraction_released_at_death = h_agent_data.fraction_released_at_death;
	d_fraction_transferred_when_ingested = h_agent_data.fraction_transferred_when_ingested;
	d_volumes = h_agent_data.volumes;

	thrust::copy_n(h_substrate_densities.get(), densities_size_bytes / sizeof(real_t), d_substrate_densities);

	// update device pointers
	{
		positions = d_positions.data().get();
		secretion_rates = d_secretion_rates.data().get();
		saturation_densities = d_saturation_densities.data().get();
		uptake_rates = d_uptake_rates.data().get();
		net_export_rates = d_net_export_rates.data().get();
		internalized_substrates = d_internalized_substrates.data().get();
		fraction_released_at_death = d_fraction_released_at_death.data().get();
		fraction_transferred_when_ingested = d_fraction_transferred_when_ingested.data().get();
		volumes = d_volumes.data().get();
	}

	residency = data_residency::DEVICE;
}
#else
void data_manager::transfer_to_host() { substrate_densities = d_substrate_densities.get(); }
void data_manager::transfer_to_device()
{
	auto& h_agent_data = retrieve_agent_data(*h_agent_container);
	positions = h_agent_data.base_data.positions.data();
	secretion_rates = h_agent_data.secretion_rates.data();
	saturation_densities = h_agent_data.saturation_densities.data();
	uptake_rates = h_agent_data.uptake_rates.data();
	net_export_rates = h_agent_data.net_export_rates.data();
	internalized_substrates = h_agent_data.internalized_substrates.data();
	fraction_released_at_death = h_agent_data.fraction_released_at_death.data();
	fraction_transferred_when_ingested = h_agent_data.fraction_transferred_when_ingested.data();
	volumes = h_agent_data.volumes.data();
}
#endif
