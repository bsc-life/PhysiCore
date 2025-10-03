#include "data_manager.h"

#include <memory>

#include "diffusion_solver.h"

using namespace physicore;
using namespace physicore::biofvm::kernels::thrust_solver;

void data_manager::initialize(microenvironment& m, diffusion_solver& d_solver)
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
void device_manager::transfer_to_host()
{
	if (residency == data_residency::HOST)
		return;

	auto& h_agent_data = retrieve_agent_data(*h_agent_container);

	h_base_agent_data.copy_to(d_base_agent_data);
	h_agent_data.copy_to(d_agent_data);
	thrust::copy_n(h_substrate_densities.get(), densities_size_bytes / sizeof(real_t), d_substrate_densities.get());

	substrate_densities = h_substrate_densities.get();

	residency = data_residency::HOST;
}

void device_manager::transfer_to_device()
{
	if (residency == data_residency::DEVICE)
		return;

	auto& h_agent_data = retrieve_agent_data(*h_agent_container);

	d_positions = h_agent_data.base_data.d_positions;
	d_secretion_rates = h_agent_data.d_secretion_rates;
	d_saturation_densities = h_agent_data.d_saturation_densities;
	d_uptake_rates = h_agent_data.d_uptake_rates;
	d_net_export_rates = h_agent_data.d_net_export_rates;
	d_internalized_substrates = h_agent_data.d_internalized_substrates;
	d_fraction_released_at_death = h_agent_data.d_fraction_released_at_death;
	d_fraction_transferred_when_ingested = h_agent_data.d_fraction_transferred_when_ingested;
	d_volumes = h_agent_data.d_volumes;
	d_substrate_densities = h_agent_data.d_substrate_densities;

	// update device pointers
	{
		positions = d_agent_data.base_data.positions.data().get();
		secretion_rates = d_agent_data.secretion_rates.data().get();
		saturation_densities = d_agent_data.saturation_densities.data().get();
		uptake_rates = d_agent_data.uptake_rates.data().get();
		net_export_rates = d_agent_data.net_export_rates.data().get();
		internalized_substrates = d_agent_data.internalized_substrates.data().get();
		fraction_released_at_death = d_agent_data.fraction_released_at_death.data().get();
		fraction_transferred_when_ingested = d_agent_data.fraction_transferred_when_ingested.data().get();
		volumes = d_agent_data.volumes.data().get();
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
