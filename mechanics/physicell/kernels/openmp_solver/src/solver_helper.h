#pragma once

#include <cmath>
#include <numbers>


namespace physicore::mechanics::physicell::kernels::openmp_solver {

constexpr void update_membrane_velocity(real_t position, real_t bounding_box, real_t sign, real_t radius,
										real_t repulsion_strength, real_t& velocity)
{
	real_t distance = std::abs(bounding_box - position);

	distance = std::max<real_t>(distance, 0.00001);

	real_t repulsion = 1 - distance / radius;
	repulsion = repulsion < 0 ? 0 : repulsion;

	repulsion *= repulsion * repulsion_strength * sign;

	velocity += repulsion * distance;
}

static constexpr void update_membrane_velocities(real_t* __restrict__ velocity, const real_t* __restrict__ position,
												 const cartesian_mesh& mesh, const real_t radius,
												 const real_t repulsion_strength)
{
	update_membrane_velocity(position[0], mesh.bounding_box_mins[0], 1, radius, repulsion_strength, velocity[0]);
	update_membrane_velocity(position[0], mesh.bounding_box_maxs[0], -1, radius, repulsion_strength, velocity[0]);
}

static constexpr void update_membrane_velocities(real_t* __restrict__ velocity, const real_t* __restrict__ position,
												 const cartesian_mesh& mesh, const real_t radius,
												 const real_t repulsion_strength)
{
	update_membrane_velocity(position[0], mesh.bounding_box_mins[0], 1, radius, repulsion_strength, velocity[0]);
	update_membrane_velocity(position[0], mesh.bounding_box_maxs[0], -1, radius, repulsion_strength, velocity[0]);
	update_membrane_velocity(position[1], mesh.bounding_box_mins[1], 1, radius, repulsion_strength, velocity[1]);
	update_membrane_velocity(position[1], mesh.bounding_box_maxs[1], -1, radius, repulsion_strength, velocity[1]);
}

static constexpr void update_membrane_velocities(real_t* __restrict__ velocity, const real_t* __restrict__ position,
												 const cartesian_mesh& mesh, const real_t radius,
												 const real_t repulsion_strength)
{
	update_membrane_velocity(position[0], mesh.bounding_box_mins[0], 1, radius, repulsion_strength, velocity[0]);
	update_membrane_velocity(position[0], mesh.bounding_box_maxs[0], -1, radius, repulsion_strength, velocity[0]);
	update_membrane_velocity(position[1], mesh.bounding_box_mins[1], 1, radius, repulsion_strength, velocity[1]);
	update_membrane_velocity(position[1], mesh.bounding_box_maxs[1], -1, radius, repulsion_strength, velocity[1]);
	update_membrane_velocity(position[2], mesh.bounding_box_mins[2], 1, radius, repulsion_strength, velocity[2]);
	update_membrane_velocity(position[2], mesh.bounding_box_maxs[2], -1, radius, repulsion_strength, velocity[2]);
}


}; // namespace physicore::mechanics::physicell::kernels::openmp_solver
