#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "microenvironment.h"
#include "microenvironment_builder.h"

using namespace physicore;
using namespace physicore::biofvm;

int main()
{
	real_t diffusion_timestep = 0.01;
	real_t simulation_time = 30.0;
	real_t output_interval = 0.1;
	index_t dims = 2;
	index_t oxygen_idx = 0;
	index_t glucose_idx = 1;

	// Create microenvironment
	std::unique_ptr<microenvironment> m;
	{
		microenvironment_builder builder;
		builder.resize(dims, { -1000, -1000, -1000 }, { 1000, 1000, 1000 }, { 20, 20, 20 });
		builder.set_time_step(diffusion_timestep);
		builder.add_density("oxygen", "mmHg", 100'000, 0.1, 38.0);
		builder.add_density("glucose", "mM", 600, 0.1, 15.0);

		oxygen_idx = builder.get_density_index("oxygen");
		glucose_idx = builder.get_density_index("glucose");

		builder.add_boundary_dirichlet_conditions(oxygen_idx, { 70.0, 60.0, 55.0 }, { 25.0, 20.0, 18.0 },
												  { true, true, true }, { true, true, true });
		builder.add_boundary_dirichlet_conditions(glucose_idx, { 9.0, 8.0, 7.5 }, { 3.0, 2.5, 2.0 },
												  { true, true, true }, { true, true, true });

		m = builder.build();
	}

	std::mt19937 rng(42);
	std::uniform_real_distribution<real_t> offset(-1.0, 1.0);

	struct agent_group
	{
		std::array<real_t, 3> center;
		real_t radius;
		std::array<real_t, 2> saturation;
		std::array<real_t, 2> uptake;
		std::array<real_t, 2> secretion;
		real_t volume;
		int count;
	};

	const std::vector<agent_group> groups = {
		{ { -600.0, 0.0, 0.0 }, 140.0, { 45.0, 6.5 }, { 500, 300 }, { 2, 1.5 }, 2200.0, 24 },
		{ { 0.0, 0.0, 0.0 }, 180.0, { 32.0, 5.0 }, { 0.09, 0.06 }, { 100, 100.5 }, 2500.0, 30 },
		{ { 600.0, 200.0, 0.0 }, 120.0, { 24.0, 4.0 }, { 0.14, 0.1 }, { 500, 500 }, 2100.0, 26 }
	};

	for (const auto& group : groups)
	{
		for (int i = 0; i < group.count; ++i)
		{
			auto* a = m->agents->create();
			auto position = a->position();
			for (std::size_t dim = 0; dim < position.size(); ++dim)
			{
				real_t displacement = offset(rng) * group.radius;
				position[dim] = group.center[dim] + displacement;
			}

			auto saturation = a->saturation_densities();
			saturation[oxygen_idx] = group.saturation[0];
			saturation[glucose_idx] = group.saturation[1];

			auto uptake = a->uptake_rates();
			uptake[oxygen_idx] = group.uptake[0];
			uptake[glucose_idx] = group.uptake[1];

			auto secretion = a->secretion_rates();
			std::fill(secretion.begin(), secretion.end(), 0.0);
			secretion[oxygen_idx] = group.secretion[0];
			secretion[glucose_idx] = group.secretion[1];

			auto net_export = a->net_export_rates();
			std::fill(net_export.begin(), net_export.end(), 0.0);

			auto internalized = a->internalized_substrates();
			std::fill(internalized.begin(), internalized.end(), 0.0);

			auto release_fraction = a->fraction_released_at_death();
			std::fill(release_fraction.begin(), release_fraction.end(), 0.0);

			auto transfer_fraction = a->fraction_transferred_when_ingested();
			std::fill(transfer_fraction.begin(), transfer_fraction.end(), 0.0);

			a->volume() = group.volume;
		}
	}

	real_t current_time = 0.0;
	real_t next_output_time = output_interval;
	std::chrono::duration<double> diffusion_runtime { 0.0 };
	std::chrono::duration<double> serialize_runtime { 0.0 };

	// m->serialize_state();

	while (current_time < simulation_time - 1e-12)
	{
		current_time += diffusion_timestep;

		auto run_start = std::chrono::steady_clock::now();
		m->run_single_timestep();
		diffusion_runtime += std::chrono::steady_clock::now() - run_start;

		if (current_time + 1e-12 >= next_output_time)
		{
			next_output_time += output_interval;

			auto serialize_start = std::chrono::steady_clock::now();
			m->serialize_state(current_time);
			serialize_runtime = std::chrono::steady_clock::now() - serialize_start;

			std::cout << "[diffuse] t=" << current_time << " diffusion runtime: " << diffusion_runtime.count() << " s"
					  << " serialization runtime: " << serialize_runtime.count() << " s" << std::endl;

			diffusion_runtime = std::chrono::duration<double> { 0.0 };
		}
	}

	return 0;
}
