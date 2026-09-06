#include <algorithm>
#include <chrono>
#include <iostream>

#include <physicell/environment.h>

#include "position_solver.h"

using namespace physicore;
using namespace physicore::mechanics::physicell;
using namespace physicore::mechanics::physicell::kernels::openmp_solver;

namespace {
void configure_agent(mechanical_agent_interface* agent, real_t radius)
{
	agent->radius() = radius;
	agent->is_movable() = 1;
	agent->is_motile() = 1;
	agent->cell_cell_adhesion_strength() = 0.4;
	agent->cell_cell_repulsion_strength() = 10.0;
	agent->cell_BM_adhesion_strength() = 4.0;
	agent->cell_BM_repulsion_strength() = 10.0;
	agent->relative_maximum_adhesion_distance() = 1.25;
	agent->maximum_number_of_attachments() = 12;
	agent->attachment_elastic_constant() = 0.01;
	agent->attachment_rate() = 0.1;
	agent->detachment_rate() = 0.1;
	agent->migration_speed() = 1.0;
	agent->migration_bias() = 0.5;
	agent->persistence_time() = 1.0;

	std::ranges::fill(agent->velocity(), 0.0);
	std::ranges::fill(agent->previous_velocity(), 0.0);
	std::ranges::fill(agent->motility_vector(), 0.0);
	std::ranges::fill(agent->migration_bias_direction(), 0.0);
	std::ranges::fill(agent->cell_adhesion_affinities(), 1.0);
}

void make_agents(environment& e, index_t count, real_t spacing)
{
	sindex_t x = 0;
	sindex_t y = 0;
	sindex_t z = 0;

	for (index_t i = 0; i < count; i++)
	{
		auto* a = e.agents->create();
		configure_agent(a, spacing * 0.5);

		auto pos = a->position();
		pos[0] = static_cast<real_t>(x);
		pos[1] = static_cast<real_t>(y);
		pos[2] = static_cast<real_t>(z);

		x += static_cast<sindex_t>(spacing);
		if (x >= e.mesh.bounding_box_maxs[0])
		{
			x -= e.mesh.bounding_box_maxs[0];
			y += static_cast<sindex_t>(spacing);
		}
		if (y >= e.mesh.bounding_box_maxs[1])
		{
			y -= e.mesh.bounding_box_maxs[1];
			z += static_cast<sindex_t>(spacing);
		}
	}
}
} // namespace

/**
 * @brief Benchmark for the physicell OpenMP position solver.
 *
 * Sets up a 3D cartesian mesh and packs it densely with mechanical agents (cells
 * touching their grid neighbors) so that neighbor search, forces, motility, basement
 * membrane interactions, spring attachments and position updates all have realistic
 * amounts of work to perform.
 *
 * The benchmark runs 100 simulation steps, measuring and printing the execution time
 * (in milliseconds) of each `position_solver` step:
 *   - update_cell_neighbors
 *   - update_cell_forces
 *   - update_motility
 *   - update_basement_membrane_interactions
 *   - update_spring_attachments
 *   - update_positions
 *
 * Timing is performed using std::chrono, and parallel execution is managed with OpenMP.
 */
int main()
{
	const real_t spacing = 20.0;
	const cartesian_mesh mesh(
		3, { 0, 0, 0 }, { 5000, 5000, 5000 },
		{ static_cast<index_t>(spacing), static_cast<index_t>(spacing), static_cast<index_t>(spacing) });

	const real_t mechanics_timestep = 0.1;

	environment e(mesh, 1, 1, mechanics_timestep);

	make_agents(e, 2'000'000, spacing);

	position_solver p_solver;

	for (index_t i = 0; i < 100; ++i)
	{
		std::size_t neighbors_duration = 0;
		std::size_t forces_duration = 0;
		std::size_t motility_duration = 0;
		std::size_t basement_membrane_duration = 0;
		std::size_t spring_attachments_duration = 0;
		std::size_t positions_duration = 0;

#pragma omp parallel private(neighbors_duration, forces_duration, motility_duration, basement_membrane_duration,       \
								 spring_attachments_duration, positions_duration)
		{
			{
				auto start = std::chrono::steady_clock::now();

				p_solver.update_cell_neighbors(e, e.mesh);

				auto end = std::chrono::steady_clock::now();

				neighbors_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			}

			{
				auto start = std::chrono::steady_clock::now();

				p_solver.update_cell_forces(e);

				auto end = std::chrono::steady_clock::now();

				forces_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			}

			{
				auto start = std::chrono::steady_clock::now();

				p_solver.update_motility(e);

				auto end = std::chrono::steady_clock::now();

				motility_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			}

			{
				auto start = std::chrono::steady_clock::now();

				p_solver.update_basement_membrane_interactions(e, e.mesh);

				auto end = std::chrono::steady_clock::now();

				basement_membrane_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			}

			{
				auto start = std::chrono::steady_clock::now();

				p_solver.update_spring_attachments(e);

				auto end = std::chrono::steady_clock::now();

				spring_attachments_duration =
					std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			}

			{
				auto start = std::chrono::steady_clock::now();

				p_solver.update_positions(e);

				auto end = std::chrono::steady_clock::now();

				positions_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			}

#pragma omp master
			std::cout << "Neighbors time: " << neighbors_duration << " ms,\t Forces time: " << forces_duration
					  << " ms,\t Motility time: " << motility_duration
					  << " ms,\t Basement membrane time: " << basement_membrane_duration
					  << " ms,\t Spring attachments time: " << spring_attachments_duration
					  << " ms,\t Positions time: " << positions_duration << " ms" << std::endl;
		}
	}
}
