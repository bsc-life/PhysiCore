#include "motility_solver.h"

#include <cmath>
#include <random>
#include <tuple>

#include <micromechanics/agent_data.h>
#include <micromechanics/environment.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void motility_solver::initialize(environment& /*e*/)
{
	if (initialized_)
		return;

	initialized_ = true;
}

void motility_solver::update_motility(environment& e)
{
	if (!e.params.enable_motility)
		return;

	auto& agents = *e.agents;
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);
	index_t const count = agents.size();
	real_t const dt = e.timestep;

	// Thread-local random generators for parallel execution
#pragma omp parallel
	{
		// Each thread gets its own random generator
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<real_t> uniform(0.0, 1.0);
		std::normal_distribution<real_t> const normal(0.0, 1.0);

#pragma omp for
		for (index_t i = 0; i < count; ++i)
		{
			if (!mech_data.is_motile[i])
				continue;

			real_t const persistence_time = mech_data.persistence_times[i];
			real_t const migration_speed = mech_data.migration_speeds[i];
			real_t const migration_bias = mech_data.migration_biases[i];

			// Check if we should update motility direction
			// Probability of update = dt / persistence_time
			if (persistence_time > 0.0 && uniform(gen) < dt / persistence_time)
			{
				// Generate random direction on unit sphere
				real_t const theta = 2.0 * M_PI * uniform(gen);
				real_t const phi = std::acos(2.0 * uniform(gen) - 1.0);

				real_t const rand_x = std::sin(phi) * std::cos(theta);
				real_t const rand_y = std::sin(phi) * std::sin(theta);
				real_t const rand_z = std::cos(phi);

				// Get bias direction
				real_t const bias_x = mech_data.migration_bias_directions[i * 3];
				real_t const bias_y = mech_data.migration_bias_directions[i * 3 + 1];
				real_t const bias_z = mech_data.migration_bias_directions[i * 3 + 2];

				// Combine random and bias: direction = (1-bias)*random + bias*bias_direction
				real_t dir_x = (1.0 - migration_bias) * rand_x + migration_bias * bias_x;
				real_t dir_y = (1.0 - migration_bias) * rand_y + migration_bias * bias_y;
				real_t dir_z = (1.0 - migration_bias) * rand_z + migration_bias * bias_z;

				// Normalize
				real_t const mag = std::sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
				if (mag > 1e-16)
				{
					dir_x /= mag;
					dir_y /= mag;
					dir_z /= mag;
				}

				// Store new motility direction
				mech_data.motility_directions[i * 3] = dir_x;
				mech_data.motility_directions[i * 3 + 1] = dir_y;
				mech_data.motility_directions[i * 3 + 2] = dir_z;
			}

			// Add motility force = speed * direction
			mech_data.forces[i * 3] += migration_speed * mech_data.motility_directions[i * 3];
			mech_data.forces[i * 3 + 1] += migration_speed * mech_data.motility_directions[i * 3 + 1];
			mech_data.forces[i * 3 + 2] += migration_speed * mech_data.motility_directions[i * 3 + 2];
		}
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
