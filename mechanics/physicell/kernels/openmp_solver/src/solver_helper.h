#pragma once

#include <cmath>
#include <numbers>

#include <common/random.h>

constexpr physicore::real_t zero_threshold = 1e-16;

namespace physicore::mechanics::physicell::kernels::openmp_solver {


constexpr void update_membrane_velocity(real_t position, real_t bounding_box, real_t sign, real_t radius,
										real_t repulsion_strength, real_t& velocity)
{
	real_t distance = bounding_box - position;
	distance = (distance < 0) ? -distance : distance;
	distance = (distance < 0.00001) ? 0.00001 : distance;

	// BITHACK VERSION: TO BE TESTED
	// distance = (distance ^ (distance >> 63)) - (distance >> 63);
	// distance += (0.00001 - distance) & -(distance < 0.00001);

	real_t repulsion = 1 - distance / radius;
	repulsion = repulsion < 0 ? 0 : repulsion;

	repulsion *= repulsion * repulsion_strength * sign;

	velocity += repulsion * distance;
}

template <index_t dims>
struct position_helper
{};

template <>
struct position_helper<1>
{
	static constexpr real_t distance(const real_t* PHYSICORE_RESTRICT lhs, const real_t* PHYSICORE_RESTRICT rhs)
	{
		return std::abs(lhs[0] - rhs[0]);
	}

	static constexpr real_t difference_and_distance(const real_t* PHYSICORE_RESTRICT lhs,
													const real_t* PHYSICORE_RESTRICT rhs,
													real_t* PHYSICORE_RESTRICT difference)
	{
		difference[0] = lhs[0] - rhs[0];

		return std::abs(difference[0]);
	}

	static constexpr void update_velocities(real_t* PHYSICORE_RESTRICT lhs, real_t* PHYSICORE_RESTRICT rhs,
											const real_t* PHYSICORE_RESTRICT difference, const real_t force)
	{
		lhs[0] += force * difference[0];
		rhs[0] -= force * difference[0];
	}

	static constexpr void update_velocity(real_t* PHYSICORE_RESTRICT velocity,
										  const real_t* PHYSICORE_RESTRICT difference, const real_t force)
	{
		velocity[0] += force * difference[0];
	}

	static void random_walk(bool, real_t* PHYSICORE_RESTRICT walk)
	{
		real_t rand = random::instance().uniform();
		walk[0] = rand < 0.5 ? -1 : 1;
	}

	static constexpr void update_motility_vector(real_t* PHYSICORE_RESTRICT motility_vector,
												 const real_t* PHYSICORE_RESTRICT walk,
												 const real_t* PHYSICORE_RESTRICT migration_bias_direction,
												 const real_t migration_bias)
	{
		motility_vector[0] = (1 - migration_bias) * walk[0] + migration_bias * migration_bias_direction[0];
	}

	static constexpr void normalize_and_scale(real_t* PHYSICORE_RESTRICT vector, real_t scale)
	{
		real_t length = std::abs(vector[0]);

		vector[0] = length > zero_threshold ? vector[0] * scale / length : 0;
	}

	static constexpr void update_membrane_velocities(real_t* PHYSICORE_RESTRICT velocity,
													 const real_t* PHYSICORE_RESTRICT position,
													 const cartesian_mesh& mesh, const real_t radius,
													 const real_t repulsion_strength)
	{
		update_membrane_velocity(position[0], mesh.bounding_box_mins[0], 1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[0], mesh.bounding_box_maxs[0], -1, radius, repulsion_strength, velocity[0]);
	}

	static constexpr void add(real_t* PHYSICORE_RESTRICT lhs, const real_t* PHYSICORE_RESTRICT rhs)
	{
		lhs[0] += rhs[0];
	}

	static constexpr void subtract(real_t* PHYSICORE_RESTRICT dst, const real_t* PHYSICORE_RESTRICT lhs,
								   const real_t* PHYSICORE_RESTRICT rhs)
	{
		dst[0] = lhs[0] - rhs[0];
	}
};

template <>
struct position_helper<2>
{
	static constexpr real_t distance(const real_t* PHYSICORE_RESTRICT lhs, const real_t* PHYSICORE_RESTRICT rhs)
	{
		return std::sqrt((lhs[0] - rhs[0]) * (lhs[0] - rhs[0]) + (lhs[1] - rhs[1]) * (lhs[1] - rhs[1]));
	}

	static constexpr real_t difference_and_distance(const real_t* PHYSICORE_RESTRICT lhs,
													const real_t* PHYSICORE_RESTRICT rhs,
													real_t* PHYSICORE_RESTRICT difference)
	{
		difference[0] = lhs[0] - rhs[0];
		difference[1] = lhs[1] - rhs[1];

		return std::sqrt(difference[0] * difference[0] + difference[1] * difference[1]);
	}

	static constexpr void update_velocities(real_t* PHYSICORE_RESTRICT lhs, real_t* PHYSICORE_RESTRICT rhs,
											const real_t* PHYSICORE_RESTRICT difference, const real_t force)
	{
		lhs[0] += force * difference[0];
		lhs[1] += force * difference[1];

		rhs[0] -= force * difference[0];
		rhs[1] -= force * difference[1];
	}

	static constexpr void update_velocity(real_t* PHYSICORE_RESTRICT velocity,
										  const real_t* PHYSICORE_RESTRICT difference, const real_t force)
	{
		velocity[0] += force * difference[0];
		velocity[1] += force * difference[1];
	}

	static void random_walk(bool, real_t* PHYSICORE_RESTRICT walk)
	{
		real_t theta = random::instance().uniform(0, 2 * std::numbers::pi_v<real_t>);
		walk[0] = std::cos(theta);
		walk[1] = std::sin(theta);
	}

	static constexpr void update_motility_vector(real_t* PHYSICORE_RESTRICT motility_vector,
												 const real_t* PHYSICORE_RESTRICT walk,
												 const real_t* PHYSICORE_RESTRICT migration_bias_direction,
												 const real_t migration_bias)
	{
		motility_vector[0] = (1 - migration_bias) * walk[0] + migration_bias * migration_bias_direction[0];
		motility_vector[1] = (1 - migration_bias) * walk[1] + migration_bias * migration_bias_direction[1];
	}

	static constexpr void normalize_and_scale(real_t* PHYSICORE_RESTRICT vector, real_t scale)
	{
		real_t length = std::sqrt(vector[0] * vector[0] + vector[1] * vector[1]);

		vector[0] = length > zero_threshold ? vector[0] * scale / length : 0;
		vector[1] = length > zero_threshold ? vector[1] * scale / length : 0;
	}

	static constexpr void update_membrane_velocities(real_t* PHYSICORE_RESTRICT velocity,
													 const real_t* PHYSICORE_RESTRICT position,
													 const cartesian_mesh& mesh, const real_t radius,
													 const real_t repulsion_strength)
	{
		update_membrane_velocity(position[0], mesh.bounding_box_mins[0], 1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[0], mesh.bounding_box_maxs[0], -1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[1], mesh.bounding_box_mins[1], 1, radius, repulsion_strength, velocity[1]);
		update_membrane_velocity(position[1], mesh.bounding_box_maxs[1], -1, radius, repulsion_strength, velocity[1]);
	}

	static constexpr void add(real_t* PHYSICORE_RESTRICT lhs, const real_t* PHYSICORE_RESTRICT rhs)
	{
		lhs[0] += rhs[0];
		lhs[1] += rhs[1];
	}

	static constexpr void subtract(real_t* PHYSICORE_RESTRICT dst, const real_t* PHYSICORE_RESTRICT lhs,
								   const real_t* PHYSICORE_RESTRICT rhs)
	{
		dst[0] = lhs[0] - rhs[0];
		dst[1] = lhs[1] - rhs[1];
	}
};

template <>
struct position_helper<3>
{
	static constexpr real_t distance(const real_t* PHYSICORE_RESTRICT lhs, const real_t* PHYSICORE_RESTRICT rhs)
	{
		return std::sqrt((lhs[0] - rhs[0]) * (lhs[0] - rhs[0]) + (lhs[1] - rhs[1]) * (lhs[1] - rhs[1])
						 + (lhs[2] - rhs[2]) * (lhs[2] - rhs[2]));
	}

	static constexpr real_t difference_and_distance(const real_t* PHYSICORE_RESTRICT lhs,
													const real_t* PHYSICORE_RESTRICT rhs,
													real_t* PHYSICORE_RESTRICT difference)
	{
		difference[0] = lhs[0] - rhs[0];
		difference[1] = lhs[1] - rhs[1];
		difference[2] = lhs[2] - rhs[2];

		return std::sqrt(difference[0] * difference[0] + difference[1] * difference[1] + difference[2] * difference[2]);
	}

	static constexpr void update_velocities(real_t* PHYSICORE_RESTRICT lhs, real_t* PHYSICORE_RESTRICT rhs,
											const real_t* PHYSICORE_RESTRICT difference, const real_t force)
	{
		lhs[0] += force * difference[0];
		lhs[1] += force * difference[1];
		lhs[2] += force * difference[2];

		rhs[0] -= force * difference[0];
		rhs[1] -= force * difference[1];
		rhs[2] -= force * difference[2];
	}

	static constexpr void update_velocity(real_t* PHYSICORE_RESTRICT velocity,
										  const real_t* PHYSICORE_RESTRICT difference, const real_t force)
	{
		velocity[0] += force * difference[0];
		velocity[1] += force * difference[1];
		velocity[2] += force * difference[2];
	}

	static void random_walk(bool restrict_to_2d, real_t* PHYSICORE_RESTRICT walk)
	{
		if (restrict_to_2d)
		{
			position_helper<2>::random_walk(true, walk);
			walk[2] = 0;
		}
		else
		{
			const real_t theta = random::instance().uniform(0, 2 * std::numbers::pi_v<real_t>);
			const real_t z = random::instance().uniform(-1, 1);
			const real_t r = std::sqrt(1 - z * z);

			walk[0] = std::cos(theta) * r;
			walk[1] = std::sin(theta) * r;
			walk[2] = z;
		}
	}

	static constexpr void update_motility_vector(real_t* PHYSICORE_RESTRICT motility_vector,
												 const real_t* PHYSICORE_RESTRICT walk,
												 const real_t* PHYSICORE_RESTRICT migration_bias_direction,
												 const real_t migration_bias)
	{
		motility_vector[0] = (1 - migration_bias) * walk[0] + migration_bias * migration_bias_direction[0];
		motility_vector[1] = (1 - migration_bias) * walk[1] + migration_bias * migration_bias_direction[1];
		motility_vector[2] = (1 - migration_bias) * walk[2] + migration_bias * migration_bias_direction[2];
	}

	static constexpr void normalize_and_scale(real_t* PHYSICORE_RESTRICT vector, real_t scale)
	{
		real_t length = std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);

		vector[0] = length > zero_threshold ? vector[0] * scale / length : 0;
		vector[1] = length > zero_threshold ? vector[1] * scale / length : 0;
		vector[2] = length > zero_threshold ? vector[2] * scale / length : 0;
	}

	static constexpr void update_membrane_velocities(real_t* PHYSICORE_RESTRICT velocity,
													 const real_t* PHYSICORE_RESTRICT position,
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

	static constexpr void add(real_t* PHYSICORE_RESTRICT lhs, const real_t* PHYSICORE_RESTRICT rhs)
	{
		lhs[0] += rhs[0];
		lhs[1] += rhs[1];
		lhs[2] += rhs[2];
	}

	static constexpr void subtract(real_t* PHYSICORE_RESTRICT dst, const real_t* PHYSICORE_RESTRICT lhs,
								   const real_t* PHYSICORE_RESTRICT rhs)
	{
		dst[0] = lhs[0] - rhs[0];
		dst[1] = lhs[1] - rhs[1];
		dst[2] = lhs[2] - rhs[2];
	}
};


}; // namespace physicore::mechanics::physicell::kernels::openmp_solver
