#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include <common/mesh.h>
#include <physicell/environment.h>

using namespace physicore;
using namespace physicore::mechanics::physicell;

namespace {

struct agent_group
{
	std::array<real_t, 3> center;
	real_t radius;
	real_t cell_radius;
	real_t adhesion_strength;
	real_t repulsion_strength;
	real_t relative_max_adhesion_distance;
	index_t count;
};

void configure_agent(mechanical_agent* agent, const agent_group& group, std::mt19937& rng)
{
	std::uniform_real_distribution<real_t> offset(-1.0, 1.0);
	auto pos = agent->position();
	for (index_t dim = 0; dim < pos.size(); ++dim)
	{
		real_t displacement = offset(rng) * group.radius;
		pos[dim] = group.center[dim] + displacement;
	}

	agent->radius() = group.cell_radius;
	agent->is_movable() = 1;
	agent->is_motile() = 0;
	agent->cell_cell_adhesion_strength() = group.adhesion_strength;
	agent->cell_cell_repulsion_strength() = group.repulsion_strength;
	agent->cell_BM_adhesion_strength() = 0.0;
	agent->cell_BM_repulsion_strength() = 30.0;
	agent->relative_maximum_adhesion_distance() = group.relative_max_adhesion_distance;
	agent->maximum_number_of_attachments() = 12;
	agent->attachment_elastic_constant() = 0.01;
	agent->attachment_rate() = 0.0;
	agent->detachment_rate() = 0.0;
	agent->migration_speed() = 0.0;
	agent->migration_bias() = 0.0;
	agent->persistence_time() = 0.0;

	auto velocity = agent->velocity();
	std::ranges::fill(velocity, 0.0);
	auto previous_velocity = agent->previous_velocity();
	std::ranges::fill(previous_velocity, 0.0);
	auto motility_vector = agent->motility_vector();
	std::ranges::fill(motility_vector, 0.0);
	auto migration_bias_direction = agent->migration_bias_direction();
	std::ranges::fill(migration_bias_direction, 0.0);
	auto orientation = agent->orientation();
	std::ranges::fill(orientation, 0.0);
	auto affinities = agent->cell_adhesion_affinities();
	std::ranges::fill(affinities, 1.0);
}

} // namespace

int main()
{
	std::unique_ptr<environment> env;
	{
		const std::filesystem::path config_file = "settings.xml";
		try
		{
			std::cout << "[mechanize] Loading configuration from: " << config_file << std::endl;
			env = environment::create_from_config(config_file);
		}
		catch (const std::exception& e)
		{
			std::cerr << "[mechanize] Error: " << e.what() << std::endl;
			return 1;
		}
	}

	std::mt19937 rng(42);
	std::uniform_real_distribution<real_t> offset(-1.0, 1.0);

	const std::vector<agent_group> groups = {
		{ .center = { -350.0, 0.0, 0.0 },
		  .radius = 90.0,
		  .cell_radius = 10.0,
		  .adhesion_strength = 0.6,
		  .repulsion_strength = 60.0,
		  .relative_max_adhesion_distance = 1.35,
		  .count = 24 },
		{ .center = { 0.0, 0.0, 0.0 },
		  .radius = 120.0,
		  .cell_radius = 11.0,
		  .adhesion_strength = 0.8,
		  .repulsion_strength = 75.0,
		  .relative_max_adhesion_distance = 1.45,
		  .count = 32 },
		{ .center = { 350.0, 150.0, 0.0 },
		  .radius = 70.0,
		  .cell_radius = 9.5,
		  .adhesion_strength = 0.7,
		  .repulsion_strength = 65.0,
		  .relative_max_adhesion_distance = 1.30,
		  .count = 20 },
		{ .center = { 0.0, -300.0, 0.0 },
		  .radius = 80.0,
		  .cell_radius = 10.5,
		  .adhesion_strength = 0.9,
		  .repulsion_strength = 80.0,
		  .relative_max_adhesion_distance = 1.40,
		  .count = 18 },
	};

	for (const auto& group : groups)
	{
		for (index_t i = 0; i < group.count; ++i)
		{
			auto* agent = env->agents->create();
			configure_agent(agent, group, rng);

			if (env->mesh.dims == 2)
			{
				agent->position()[0] += offset(rng) * group.radius * 0.25;
				agent->position()[1] += offset(rng) * group.radius * 0.25;
			}
			else
			{
				agent->position()[0] += offset(rng) * group.radius * 0.25;
				agent->position()[1] += offset(rng) * group.radius * 0.25;
				agent->position()[2] += offset(rng) * group.radius * 0.25;
			}
		}
	}

	std::cout << "[mechanize] Created " << env->agents->size() << " cells across " << groups.size()
			  << " dense groups to stress cell-cell and wall interactions." << std::endl;

	const real_t output_interval = 0.1;
	real_t current_time = 0.0;
	real_t next_output_time = output_interval;
	std::chrono::duration<double> mechanics_runtime { 0.0 };
	std::chrono::duration<double> serialize_runtime { 0.0 };

	env->serialize_state(current_time);

	std::cout << "\n[mechanize] Running simulation for " << env->simulation_time << " time units..." << std::endl;

	while (current_time < env->simulation_time - 1e-12)
	{
		current_time += env->mechanics_timestep;

		auto run_start = std::chrono::steady_clock::now();
		env->run_single_timestep();
		mechanics_runtime += std::chrono::steady_clock::now() - run_start;

		if (current_time + 1e-12 >= next_output_time)
		{
			next_output_time += output_interval;

			auto serialize_start = std::chrono::steady_clock::now();
			env->serialize_state(current_time);
			serialize_runtime = std::chrono::steady_clock::now() - serialize_start;

			std::cout << "[mechanize] t=" << current_time << " mechanics runtime: " << mechanics_runtime.count() << " s"
					  << " serialization runtime: " << serialize_runtime.count() << " s" << std::endl;

			mechanics_runtime = std::chrono::duration<double> { 0.0 };
		}
	}

	std::cout << "\n[mechanize] Simulation completed successfully!" << std::endl;
	return 0;
}
