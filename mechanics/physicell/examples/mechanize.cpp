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
#include <physicell/solver_registry.h>
#include <physicell/vtk_agents_serializer.h>

#include "../src/config_reader.h"

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

std::array<sindex_t, 3> to_min_bounds(const domain_config& domain)
{
	return { static_cast<sindex_t>(std::lround(domain.x_min)), static_cast<sindex_t>(std::lround(domain.y_min)),
			 static_cast<sindex_t>(std::lround(domain.z_min)) };
}

std::array<sindex_t, 3> to_max_bounds(const domain_config& domain)
{
	return { static_cast<sindex_t>(std::lround(domain.x_max)), static_cast<sindex_t>(std::lround(domain.y_max)),
			 static_cast<sindex_t>(std::lround(domain.z_max)) };
}

std::array<index_t, 3> to_voxel_shape(const domain_config& domain)
{
	return { static_cast<index_t>(std::max<real_t>(1.0, domain.dx)),
			 static_cast<index_t>(std::max<real_t>(1.0, domain.dy)),
			 static_cast<index_t>(std::max<real_t>(1.0, domain.dz)) };
}

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
	auto resolve_config = []() -> std::filesystem::path {
		const std::array<std::filesystem::path, 4> candidates = {
			std::filesystem::path { "settings.xml" },
			std::filesystem::path { "mechanics/physicell/examples/settings.xml" },
			std::filesystem::path { "/workspaces/PhysiCore/mechanics/physicell/examples/settings.xml" },
			std::filesystem::current_path() / "mechanics/physicell/examples/settings.xml"
		};

		for (const auto& candidate : candidates)
		{
			if (std::filesystem::exists(candidate))
			{
				return candidate;
			}
		}

		std::filesystem::path executable_dir = std::filesystem::path { "/proc/self/exe" }.parent_path();
		if (std::filesystem::exists(executable_dir / "settings.xml"))
		{
			return executable_dir / "settings.xml";
		}

		return "settings.xml";
	};

	const std::filesystem::path config_file = resolve_config();
	mechanics_config config;

	try
	{
		std::cout << "[mechanize] Loading configuration from: " << config_file << std::endl;
		config = parse_simulation_parameters(config_file);
	}
	catch (const std::exception& e)
	{
		std::cerr << "[mechanize] Error: " << e.what() << std::endl;
		return 1;
	}

	const index_t dims = config.is_2D ? 2 : 3;
	environment env(config.overall.dt_mechanics, dims, static_cast<index_t>(config.cell_types.size()), 0);
	env.set_mesh(cartesian_mesh { dims, to_min_bounds(config.domain), to_max_bounds(config.domain),
								  to_voxel_shape(config.domain) });

	auto solver = solver_registry::instance().get("openmp_solver");
	if (!solver)
	{
		std::cerr << "[mechanize] Error: openmp_solver not registered" << std::endl;
		return 2;
	}
	env.solver = std::move(solver);
	env.solver->initialize(env);

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
			auto* agent = env.agents->create();
			configure_agent(agent, group, rng);

			if (dims == 2)
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

	std::cout << "[mechanize] Created " << env.agents->size() << " cells across " << groups.size()
			  << " dense groups to stress cell-cell and wall interactions." << std::endl;

	const std::vector<std::string> substrate_names = { "oxygen", "necrotic debris", "apoptotic debris" };
	const std::vector<std::string> cell_type_names = { "malignant epithelial" };
	env.serializer = std::make_unique<vtk_agents_serializer>(
		"mechanics_vtk_output", *static_cast<mechanical_agent_container_interface*>(env.agents.get()), substrate_names,
		cell_type_names);

	const real_t output_interval = 0.1;
	real_t current_time = 0.0;
	real_t next_output_time = output_interval;
	std::chrono::duration<double> mechanics_runtime { 0.0 };
	std::chrono::duration<double> serialize_runtime { 0.0 };

	env.serialize_state(current_time);

	std::cout << "\n[mechanize] Running simulation for " << config.overall.max_time << " time units..."
			  << std::endl;

	while (current_time < config.overall.max_time - 1e-12)
	{
		current_time += config.overall.dt_mechanics;

		auto run_start = std::chrono::steady_clock::now();
		env.run_single_timestep();
		mechanics_runtime += std::chrono::steady_clock::now() - run_start;

		if (current_time + 1e-12 >= next_output_time)
		{
			next_output_time += output_interval;

			auto serialize_start = std::chrono::steady_clock::now();
			env.serialize_state(current_time);
			serialize_runtime = std::chrono::steady_clock::now() - serialize_start;

			std::cout << "[mechanize] t=" << current_time << " mechanics runtime: " << mechanics_runtime.count()
				  << " s"
				  << " serialization runtime: " << serialize_runtime.count() << " s" << std::endl;

			mechanics_runtime = std::chrono::duration<double> { 0.0 };
		}
	}

	std::cout << "\n[mechanize] Simulation completed successfully!" << std::endl;
	return 0;
}
