#include <chrono>
#include <iostream>

#include <biofvm/bulk_functor.h>
#include <biofvm/microenvironment.h>

#include "bulk_solver.h"
#include "cell_solver.h"
#include "diffusion_solver.h"
#include "dirichlet_solver.h"

using namespace physicore;
using namespace physicore::reactions_diffusion::biofvm;
using namespace physicore::reactions_diffusion::biofvm::kernels::openmp_solver;

namespace {
void make_agents(microenvironment& m, index_t count, bool conflict)
{
	sindex_t x = 0;
	sindex_t y = 0;
	sindex_t z = 0;

	for (index_t i = 0; i < count; i++)
	{
		auto* a = m.agents->create();
		a->position()[0] = static_cast<real_t>(x);
		a->position()[1] = static_cast<real_t>(y);
		a->position()[2] = static_cast<real_t>(z);

		x += 20;
		if (x >= m.mesh.bounding_box_maxs[0])
		{
			x -= m.mesh.bounding_box_maxs[0];
			y += 20;
		}
		if (y >= m.mesh.bounding_box_maxs[1])
		{
			y -= m.mesh.bounding_box_maxs[1];
			z += 20;
		}
	}

	if (conflict)
	{
		auto* a = m.agents->create();
		a->position()[0] = 0;
		a->position()[1] = 0;
		a->position()[2] = 0;
	}
}
} // namespace

/**
 * @brief Benchmark for reaction-diffusion solvers in a microenvironment.
 *
 * This main function sets up a 3D cartesian mesh and initializes a microenvironment
 * with multiple substrates. It configures initial conditions, diffusion coefficients,
 * decay rates, and boundary conditions for the simulation. Bulk supply, uptake, and
 * target density functions are provided as simple lambdas for benchmarking purposes.
 *
 * The function creates a large number of agents and initializes the bulk, cell,
 * diffusion, and Dirichlet solvers. It then runs 100 simulation steps, measuring
 * and printing the execution time (in milliseconds) for each solver step:
 *   - Diffusion solver
 *   - Secretion and uptake solver
 *   - Bulk solver
 *   - Dirichlet boundary condition solver
 *
 * Timing is performed using std::chrono, and parallel execution is managed with OpenMP.
 *
 * This benchmark is intended to evaluate the performance of the solvers under
 * realistic simulation conditions.
 */
int main()
{
	const cartesian_mesh mesh(3, { 0, 0, 0 }, { 5000, 5000, 5000 }, { 20, 20, 20 });

	const real_t diffusion_timestep = 0.01;
	const index_t substrates_count = 4;

	microenvironment m(mesh, substrates_count, diffusion_timestep);
	m.compute_internalized_substrates = true;

	// --- diffusion parameters -------------------------------
	m.initial_conditions = std::make_unique<real_t[]>(substrates_count);
	m.diffusion_coefficients = std::make_unique<real_t[]>(substrates_count);
	m.decay_rates = std::make_unique<real_t[]>(substrates_count);
	for (index_t i = 0; i < substrates_count; i++)
	{
		m.initial_conditions[i] = 1000;
		m.diffusion_coefficients[i] = 10000;
		m.decay_rates[i] = 0.5;
	}

	// --- bulk functors -------------------------------------------------
	struct benchmark_bulk_functor : bulk_functor
	{
		real_t supply_rates(index_t s, index_t /*x*/, index_t /*y*/, index_t /*z*/) override
		{
			return 0.01 * static_cast<real_t>(s + 1);
		}

		real_t supply_target_densities(index_t s, index_t /*x*/, index_t /*y*/, index_t /*z*/) override
		{
			return 0.01 * static_cast<real_t>(s + 1);
		}

		real_t uptake_rates(index_t s, index_t /*x*/, index_t /*y*/, index_t /*z*/) override
		{
			return 0.01 * static_cast<real_t>(s + 1);
		}
	};
	m.bulk_fnc = std::make_unique<benchmark_bulk_functor>();

	// --- dirichlet / boundary conditions -------------------------------
	for (index_t dim = 0; dim < m.mesh.dims; ++dim)
	{
		m.dirichlet_min_boundary_values[dim] = std::make_unique<real_t[]>(m.substrates_count);
		m.dirichlet_max_boundary_values[dim] = std::make_unique<real_t[]>(m.substrates_count);
		m.dirichlet_min_boundary_conditions[dim] = std::make_unique<bool[]>(m.substrates_count);
		m.dirichlet_max_boundary_conditions[dim] = std::make_unique<bool[]>(m.substrates_count);

		for (index_t s = 0; s < m.substrates_count; ++s)
		{
			m.dirichlet_min_boundary_values[dim][s] = 0.0;
			m.dirichlet_max_boundary_values[dim][s] = 0.5 * static_cast<real_t>(dim + 1);

			m.dirichlet_min_boundary_conditions[dim][s] = true;
			m.dirichlet_max_boundary_conditions[dim][s] = true;
		}
	}

	make_agents(m, 2'000'000, true);


	bulk_solver b_solver;
	cell_solver c_solver;
	diffusion_solver d_solver;
	// dirichlet_solver has only static API in this example; no instance needed

	b_solver.initialize(m);
	c_solver.initialize(m);
	d_solver.prepare(m, 1);
	d_solver.initialize();


	for (index_t i = 0; i < 100; ++i)
	{
		std::size_t diffusion_duration = 0;
		std::size_t secretion_duration = 0;
		std::size_t bulk_duration = 0;
		std::size_t dirichlet_duration = 0;

#pragma omp parallel private(diffusion_duration, secretion_duration, bulk_duration, dirichlet_duration)
		{
			{
				auto start = std::chrono::steady_clock::now();

				d_solver.solve();

				auto end = std::chrono::steady_clock::now();

				diffusion_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			}

			{
				auto start = std::chrono::steady_clock::now();

				c_solver.simulate_secretion_and_uptake(m, d_solver, i % 10 == 0);

				auto end = std::chrono::steady_clock::now();

				secretion_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			}

			// bulk solver benchmark
			{
				auto start = std::chrono::steady_clock::now();

				b_solver.solve(m, d_solver);

				auto end = std::chrono::steady_clock::now();

				bulk_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			}

			// dirichlet solver benchmark
			{
				auto start = std::chrono::steady_clock::now();

				dirichlet_solver::solve(m, d_solver);

				auto end = std::chrono::steady_clock::now();

				dirichlet_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			}

#pragma omp master
			std::cout << "Diffusion time: " << diffusion_duration << " ms,\t Secretion time: " << secretion_duration
					  << " ms,\t Bulk time: " << bulk_duration << " ms,\t Dirichlet time: " << dirichlet_duration
					  << " ms" << std::endl;
		}
	}
}
