#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "microenvironment.h"

using namespace physicore;
using namespace physicore::biofvm;

int main()
{
	real_t output_interval = 0.1;
	index_t oxygen_idx = 0;
	index_t glucose_idx = 1;

	// Create microenvironment
	std::unique_ptr<microenvironment> m;
	{
		std::filesystem::path config_file = "settings.xml";

		try
		{
			std::cout << "[diffuse] Loading configuration from: " << config_file << std::endl;
			m = microenvironment::create_from_config(config_file);
		}
		catch (const std::exception& e)
		{
			std::cerr << "[diffuse] Error: " << e.what() << std::endl;
			return 1;
		}
	}

	m->print_info(std::cout);

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

	m->solver->initialize(*m);
	m->serialize_state(current_time);

	std::cout << "\n[diffuse] Running simulation for " << m->simulation_time << " time units..." << std::endl;

	while (current_time < m->simulation_time - 1e-12)
	{
		current_time += m->diffusion_timestep;

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

	std::cout << "\n[diffuse] Simulation completed successfully!" << std::endl;

	return 0;
}
