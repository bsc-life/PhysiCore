#include "../include/physicell/openmp_solver/position_solver.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <common/cartesian_mesh.h>
#include <common/types.h>
#include <physicell/mechanical_agent_container.h>

#include "common_solver.h"
#include "random.h"
#include "solver_helper.h"

using namespace physicore::mechanics::physicell;

namespace physicore::mechanics::physicell::kernels::openmp_solver {
constexpr real_t simple_pressure_coefficient = 36.64504274775163; // 1 / (12 * (1 - sqrt(pi/(2*sqrt(3))))^2)

namespace {

void clear_simple_pressure(real_t* PHYSICORE_RESTRICT simple_pressure, index_t count)
{
#pragma omp for
	for (index_t i = 0; i < count; i++)
	{
		simple_pressure[i] = 0;
	}
}

template <index_t dims>
void solve_pair(index_t lhs, index_t rhs, index_t cell_defs_count, real_t* PHYSICORE_RESTRICT velocity,
				real_t* PHYSICORE_RESTRICT simple_pressure, const real_t* PHYSICORE_RESTRICT position,
				const real_t* PHYSICORE_RESTRICT radius, const real_t* PHYSICORE_RESTRICT cell_cell_repulsion_strength,
				const real_t* PHYSICORE_RESTRICT cell_cell_adhesion_strength,
				const real_t* PHYSICORE_RESTRICT relative_maximum_adhesion_distance,
				const real_t* PHYSICORE_RESTRICT cell_adhesion_affinity,
				const index_t* PHYSICORE_RESTRICT cell_definition_index)
{
	std::array<real_t, dims> position_difference {};

	const real_t distance =
		std::max<real_t>(position_helper<dims>::difference_and_distance(position + lhs * dims, position + rhs * dims,
																		position_difference.data()),
						 0.00001);

	// compute repulsion
	real_t repulsion {};
	{
		const real_t repulsive_distance = radius[lhs] + radius[rhs];

		repulsion = 1 - distance / repulsive_distance;

		repulsion = repulsion < 0 ? 0 : repulsion;

		repulsion *= repulsion;

		// update simple pressure
		simple_pressure[lhs] += repulsion * simple_pressure_coefficient;
		simple_pressure[rhs] += repulsion * simple_pressure_coefficient;

		repulsion *= std::sqrt(cell_cell_repulsion_strength[lhs] * cell_cell_repulsion_strength[rhs]);
	}

	// compute adhesion
	real_t adhesion {};
	{
		const real_t adhesion_distance = relative_maximum_adhesion_distance[lhs] * radius[lhs]
										 + relative_maximum_adhesion_distance[rhs] * radius[rhs];

		adhesion = 1 - distance / adhesion_distance;

		adhesion = adhesion < 0 ? 0 : adhesion;

		adhesion *= adhesion;

		const index_t lhs_cell_def_index = cell_definition_index[lhs];
		const index_t rhs_cell_def_index = cell_definition_index[rhs];

		adhesion *= std::sqrt(cell_cell_adhesion_strength[lhs] * cell_cell_adhesion_strength[rhs]
							  * cell_adhesion_affinity[lhs * cell_defs_count + rhs_cell_def_index]
							  * cell_adhesion_affinity[rhs * cell_defs_count + lhs_cell_def_index]);
	}

	real_t const force = (repulsion - adhesion) / distance;

	position_helper<dims>::update_velocity(velocity + lhs * dims, position_difference.data(), force);
}

template <index_t dims>
void update_cell_forces_single(index_t i, index_t cell_def_count, real_t* PHYSICORE_RESTRICT velocity,
							   real_t* PHYSICORE_RESTRICT simple_pressure, const real_t* PHYSICORE_RESTRICT position,
							   const real_t* PHYSICORE_RESTRICT radius,
							   const real_t* PHYSICORE_RESTRICT cell_cell_repulsion_strength,
							   const real_t* PHYSICORE_RESTRICT cell_cell_adhesion_strength,
							   const real_t* PHYSICORE_RESTRICT relative_maximum_adhesion_distance,
							   const index_t* PHYSICORE_RESTRICT cell_definition_index,
							   const real_t* PHYSICORE_RESTRICT cell_adhesion_affinities,
							   std::vector<index_t>* PHYSICORE_RESTRICT neighbors)
{
	for (const index_t j : neighbors[i])
	{
		solve_pair<dims>(i, j, cell_def_count, velocity, simple_pressure, position, radius,
						 cell_cell_repulsion_strength, cell_cell_adhesion_strength, relative_maximum_adhesion_distance,
						 cell_adhesion_affinities, cell_definition_index);
	}
}

template <index_t dims>
void update_cell_forces_internal(index_t agents_count, index_t cell_def_count, real_t* PHYSICORE_RESTRICT velocity,
								 real_t* PHYSICORE_RESTRICT simple_pressure, const real_t* PHYSICORE_RESTRICT position,
								 const real_t* PHYSICORE_RESTRICT radius,
								 const real_t* PHYSICORE_RESTRICT cell_cell_repulsion_strength,
								 const real_t* PHYSICORE_RESTRICT cell_cell_adhesion_strength,
								 const real_t* PHYSICORE_RESTRICT relative_maximum_adhesion_distance,
								 const index_t* PHYSICORE_RESTRICT cell_definition_index,
								 const real_t* PHYSICORE_RESTRICT cell_adhesion_affinities,
								 const std::uint8_t* PHYSICORE_RESTRICT is_movable,
								 std::vector<index_t>* PHYSICORE_RESTRICT neighbors)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		if (is_movable[i] == 0)
			continue;

		update_cell_forces_single<dims>(i, cell_def_count, velocity, simple_pressure, position, radius,
										cell_cell_repulsion_strength, cell_cell_adhesion_strength,
										relative_maximum_adhesion_distance, cell_definition_index,
										cell_adhesion_affinities, neighbors);
	}
}

template <index_t dims>
void update_cell_neighbors_single(environment& e, index_t i, const real_t* PHYSICORE_RESTRICT position,
								  const real_t* PHYSICORE_RESTRICT radius,
								  const real_t* PHYSICORE_RESTRICT relative_maximum_adhesion_distance,
								  std::vector<index_t>* PHYSICORE_RESTRICT neighbors, const cartesian_mesh& mesh,
								  const std::vector<std::vector<index_t>>& cells_in_voxels)
{
	(void)e;
	common_solver::for_each_in_mech_neighborhood(
		mesh, cells_in_voxels, common_solver::get_mesh_position(position + dims * i, mesh), i, [=](index_t j) {
			const real_t adhesion_distance =
				relative_maximum_adhesion_distance[i] * radius[i] + relative_maximum_adhesion_distance[j] * radius[j];

			const real_t distance = position_helper<dims>::distance(position + i * dims, position + j * dims);

			if (distance <= adhesion_distance)
			{
				neighbors[i].push_back(j);
			}
		});
}

template <index_t dims>
void update_cell_neighbors_internal(environment& e, index_t agents_count, const real_t* PHYSICORE_RESTRICT position,
									const real_t* PHYSICORE_RESTRICT radius,
									const real_t* PHYSICORE_RESTRICT relative_maximum_adhesion_distance,
									const std::uint8_t* PHYSICORE_RESTRICT is_movable,
									std::vector<index_t>* PHYSICORE_RESTRICT neighbors, const cartesian_mesh& mesh,
									const std::vector<std::vector<index_t>>& cells_in_voxels)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		if (is_movable[i] == 0)
			continue;

		update_cell_neighbors_single<dims>(e, i, position, radius, relative_maximum_adhesion_distance, neighbors, mesh,
										   cells_in_voxels);
	}
}

template <index_t dims>
void update_motility_single(
	index_t i, real_t time_step, real_t* PHYSICORE_RESTRICT motility_vector, real_t* PHYSICORE_RESTRICT velocity,
	const real_t* PHYSICORE_RESTRICT persistence_time, const real_t* PHYSICORE_RESTRICT migration_bias,
	real_t* PHYSICORE_RESTRICT migration_bias_direction, const std::uint8_t* PHYSICORE_RESTRICT restrict_to_2d,
	const std::uint8_t* PHYSICORE_RESTRICT is_motile, const real_t* PHYSICORE_RESTRICT migration_speed,
	const motility_properties::direction_update_func* PHYSICORE_RESTRICT update_migration_bias_direction_f,
	const index_t* PHYSICORE_RESTRICT cell_definition_index)
{
	if (is_motile[i] == 0)
		return;

	if (physicore::random::uniform() < time_step / persistence_time[i])
	{
		std::array<real_t, dims> random_walk {};

		position_helper<dims>::random_walk(restrict_to_2d, random_walk.data());

		if (update_migration_bias_direction_f != nullptr && update_migration_bias_direction_f[i])
		{
			update_migration_bias_direction_f[i](cell_definition_index[i]);
		}

		position_helper<dims>::update_motility_vector(motility_vector + i * dims, random_walk.data(),
													  migration_bias_direction + i * dims, migration_bias[i]);

		position_helper<dims>::normalize_and_scale(motility_vector + i * dims, migration_speed[i]);
	}

	position_helper<dims>::add(velocity + i * dims, motility_vector + i * dims);
}

template <index_t dims>
void update_motility_internal(
	index_t agents_count, real_t time_step, real_t* PHYSICORE_RESTRICT motility_vector,
	real_t* PHYSICORE_RESTRICT velocity, const real_t* PHYSICORE_RESTRICT persistence_time,
	const real_t* PHYSICORE_RESTRICT migration_bias, real_t* PHYSICORE_RESTRICT migration_bias_direction,
	const std::uint8_t* PHYSICORE_RESTRICT restrict_to_2d, const std::uint8_t* PHYSICORE_RESTRICT is_motile,
	const real_t* PHYSICORE_RESTRICT migration_speed,
	const motility_properties::direction_update_func* PHYSICORE_RESTRICT update_migration_bias_direction_f,
	const index_t* PHYSICORE_RESTRICT cell_definition_index)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		update_motility_single<dims>(i, time_step, motility_vector, velocity, persistence_time, migration_bias,
									 migration_bias_direction, restrict_to_2d, is_motile, migration_speed,
									 update_migration_bias_direction_f, cell_definition_index);
	}
}

template <index_t dims>
void update_basement_membrane_interactions_single(index_t i, real_t* PHYSICORE_RESTRICT velocity,
												  const real_t* PHYSICORE_RESTRICT position,
												  const real_t* PHYSICORE_RESTRICT radius,
												  const real_t* PHYSICORE_RESTRICT cell_BM_repulsion_strength,
												  const cartesian_mesh& mesh)
{
	position_helper<dims>::update_membrane_velocities(velocity + i * dims, position + i * dims, mesh, radius[i],
													  cell_BM_repulsion_strength[i]);
}

template <index_t dims>
void update_basement_membrane_interactions_internal(index_t agents_count, real_t* PHYSICORE_RESTRICT velocity,
													const real_t* PHYSICORE_RESTRICT position,
													const real_t* PHYSICORE_RESTRICT radius,
													const real_t* PHYSICORE_RESTRICT cell_BM_repulsion_strength,
													const std::uint8_t* PHYSICORE_RESTRICT is_movable,
													const cartesian_mesh& mesh)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		if (is_movable[i] == 0)
			continue;

		update_basement_membrane_interactions_single<dims>(i, velocity, position, radius, cell_BM_repulsion_strength,
														   mesh);
	}
}

template <index_t dims>
void spring_contract_function(index_t agents_count, index_t cell_defs_count, real_t* PHYSICORE_RESTRICT velocity,
							  const index_t* PHYSICORE_RESTRICT cell_definition_index,
							  const real_t* PHYSICORE_RESTRICT attachment_elastic_constant,
							  const real_t* PHYSICORE_RESTRICT cell_adhesion_affinity,
							  const real_t* PHYSICORE_RESTRICT position,
							  const std::uint8_t* PHYSICORE_RESTRICT is_movable,
							  std::vector<index_t>* PHYSICORE_RESTRICT springs)
{
#pragma omp for
	for (index_t this_cell_index = 0; this_cell_index < agents_count; this_cell_index++)
	{
		if (is_movable[this_cell_index] == 0)
			continue;

		for (std::size_t j = 0; j < springs[this_cell_index].size(); j++)
		{
			const index_t other_cell_index = springs[this_cell_index][j];

			const index_t this_cell_def_index = cell_definition_index[this_cell_index];
			const index_t other_cell_def_index = cell_definition_index[other_cell_index];

			const real_t adhesion =
				sqrt(attachment_elastic_constant[this_cell_index] * attachment_elastic_constant[other_cell_index]
					 * cell_adhesion_affinity[this_cell_index * cell_defs_count + other_cell_def_index]
					 * cell_adhesion_affinity[other_cell_index * cell_defs_count + this_cell_def_index]);

			std::array<real_t, dims> difference {};

			position_helper<dims>::subtract(difference.data(), position + other_cell_index * dims,
											position + this_cell_index * dims);

			position_helper<dims>::update_velocity(velocity + this_cell_index * dims, difference.data(), adhesion);
		}
	}
}

constexpr index_t erased_spring = -1;

void update_spring_attachments_internal(index_t agents_count, real_t time_step, index_t cell_defs_count,
										const real_t* PHYSICORE_RESTRICT detachment_rate,
										const real_t* PHYSICORE_RESTRICT attachment_rate,
										const real_t* PHYSICORE_RESTRICT cell_adhesion_affinities,
										const index_t* PHYSICORE_RESTRICT maximum_number_of_attachments,
										const index_t* PHYSICORE_RESTRICT cell_definition_index,
										const std::vector<index_t>* PHYSICORE_RESTRICT neighbors,
										std::vector<index_t>* PHYSICORE_RESTRICT springs)
{
// mark springs for detachment
#pragma omp for
	for (index_t this_cell_index = 0; this_cell_index < agents_count; this_cell_index++)
	{
		for (index_t j = 0; j < (index_t)springs[this_cell_index].size(); j++)
		{
			if (physicore::random::uniform() <= detachment_rate[this_cell_index] * time_step)
			{
#pragma omp critical
				{
					const index_t other_cell_index = springs[this_cell_index][j];

					if (other_cell_index != erased_spring)
					{
						springs[this_cell_index][j] = erased_spring;

						*std::ranges::find(springs[other_cell_index], this_cell_index) = erased_spring;
					}
				}
			}
		}
	}

// remove marked springs
#pragma omp for
	for (index_t this_cell_index = 0; this_cell_index < agents_count; this_cell_index++)
	{
		auto removed = std::ranges::remove(springs[this_cell_index], erased_spring);

		springs[this_cell_index].erase(removed.begin(), removed.end());
	}

	// attach cells to springs

#pragma omp for
	for (index_t this_cell_index = 0; this_cell_index < agents_count; this_cell_index++)
	{
		for (std::size_t j = 0; j < neighbors[this_cell_index].size(); j++)
		{
			const index_t other_cell_index = neighbors[this_cell_index][j];

			if (other_cell_index < this_cell_index)
				continue;

			const real_t affinity_l =
				cell_adhesion_affinities[this_cell_index * cell_defs_count + cell_definition_index[other_cell_index]];

			const real_t attachment_prob_l = attachment_rate[this_cell_index] * time_step * affinity_l;

			const real_t affinity_r =
				cell_adhesion_affinities[other_cell_index * cell_defs_count + cell_definition_index[this_cell_index]];

			const real_t attachment_prob_r = attachment_rate[other_cell_index] * time_step * affinity_r;

			if (physicore::random::uniform() <= attachment_prob_l || physicore::random::uniform() <= attachment_prob_r)
			{
#pragma omp critical
				{
					if ((index_t)springs[this_cell_index].size() < maximum_number_of_attachments[this_cell_index]
						&& (index_t)springs[other_cell_index].size() < maximum_number_of_attachments[other_cell_index])
					{
						springs[this_cell_index].push_back(other_cell_index);
						springs[other_cell_index].push_back(this_cell_index);
					}
				}
			}
		}
	}
}

template <index_t dims>
void update_positions_internal(index_t agents_count, real_t time_step, real_t* PHYSICORE_RESTRICT position,
							   real_t* PHYSICORE_RESTRICT velocity, real_t* PHYSICORE_RESTRICT previous_velocity,
							   const std::uint8_t* PHYSICORE_RESTRICT is_movable)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		if (!is_movable[i])
			continue;

		const real_t factor = time_step * 1.5;
		const real_t previous_factor = time_step * -0.5;

		for (index_t d = 0; d < dims; d++)
		{
			position[i * dims + d] +=
				velocity[i * dims + d] * factor + previous_velocity[i * dims + d] * previous_factor;

			previous_velocity[i * dims + d] = velocity[i * dims + d];
			velocity[i * dims + d] = 0;
		}
	}
}

} // namespace

void position_solver::update_cell_forces(environment& e)
{
	auto& data = retrieve_agent_data(*e.agents);
	const index_t dims = data.base_data.dims;

	clear_simple_pressure(data.state_data.simple_pressure.data(), data.agents_count);

	if (dims == 1)
		update_cell_forces_internal<1>(data.agents_count, data.agent_types_count, data.velocity.data(),
									   data.state_data.simple_pressure.data(), data.base_data.positions.data(),
									   data.radius.data(), data.mechanics_data.cell_cell_repulsion_strength.data(),
									   data.mechanics_data.cell_cell_adhesion_strength.data(),
									   data.mechanics_data.relative_maximum_adhesion_distance.data(),
									   data.state_data.agent_type_index.data(),
									   data.mechanics_data.cell_adhesion_affinities.data(),
									   data.state_data.is_movable.data(), data.state_data.neighbors.data());
	else if (dims == 2)
		update_cell_forces_internal<2>(data.agents_count, data.agent_types_count, data.velocity.data(),
									   data.state_data.simple_pressure.data(), data.base_data.positions.data(),
									   data.radius.data(), data.mechanics_data.cell_cell_repulsion_strength.data(),
									   data.mechanics_data.cell_cell_adhesion_strength.data(),
									   data.mechanics_data.relative_maximum_adhesion_distance.data(),
									   data.state_data.agent_type_index.data(),
									   data.mechanics_data.cell_adhesion_affinities.data(),
									   data.state_data.is_movable.data(), data.state_data.neighbors.data());
	else if (dims == 3)
		update_cell_forces_internal<3>(data.agents_count, data.agent_types_count, data.velocity.data(),
									   data.state_data.simple_pressure.data(), data.base_data.positions.data(),
									   data.radius.data(), data.mechanics_data.cell_cell_repulsion_strength.data(),
									   data.mechanics_data.cell_cell_adhesion_strength.data(),
									   data.mechanics_data.relative_maximum_adhesion_distance.data(),
									   data.state_data.agent_type_index.data(),
									   data.mechanics_data.cell_adhesion_affinities.data(),
									   data.state_data.is_movable.data(), data.state_data.neighbors.data());
}

void position_solver::update_cell_neighbors(environment& e, const cartesian_mesh& mesh)
{
	auto& data = retrieve_agent_data(*e.agents);
	const index_t dims = data.base_data.dims;

	std::vector<std::vector<index_t>> cells_in_voxels(mesh.voxel_count());
	for (index_t i = 0; i < data.agents_count; i++)
	{
		const auto voxel_pos = common_solver::get_mesh_position(data.base_data.positions.data() + dims * i, mesh);
		const auto voxel_idx = common_solver::get_mesh_index(voxel_pos, mesh);
		cells_in_voxels[voxel_idx].push_back(i);
	}

	// clear neighbors
#pragma omp for
	for (index_t i = 0; i < data.agents_count; i++)
		data.state_data.neighbors[i].clear();

	if (dims == 1)
		update_cell_neighbors_internal<1>(e, data.agents_count, data.base_data.positions.data(), data.radius.data(),
										  data.mechanics_data.relative_maximum_adhesion_distance.data(),
										  data.state_data.is_movable.data(), data.state_data.neighbors.data(), mesh,
										  cells_in_voxels);
	else if (dims == 2)
		update_cell_neighbors_internal<2>(e, data.agents_count, data.base_data.positions.data(), data.radius.data(),
										  data.mechanics_data.relative_maximum_adhesion_distance.data(),
										  data.state_data.is_movable.data(), data.state_data.neighbors.data(), mesh,
										  cells_in_voxels);
	else if (dims == 3)
		update_cell_neighbors_internal<3>(e, data.agents_count, data.base_data.positions.data(), data.radius.data(),
										  data.mechanics_data.relative_maximum_adhesion_distance.data(),
										  data.state_data.is_movable.data(), data.state_data.neighbors.data(), mesh,
										  cells_in_voxels);
}

void position_solver::update_motility(environment& e)
{
	auto& data = retrieve_agent_data(*e.agents);

	if (data.base_data.dims == 1)
		update_motility_internal<1>(
			data.agents_count, e.mechanics_timestep, data.motility_data.motility_vector.data(), data.velocity.data(),
			data.motility_data.persistence_time.data(), data.motility_data.migration_bias.data(),
			data.motility_data.migration_bias_direction.data(), data.motility_data.restrict_to_2d.data(),
			data.motility_data.is_motile.data(), data.motility_data.migration_speed.data(),
			data.motility_data.direction_update_funcs.data(), data.state_data.agent_type_index.data());
	else if (data.base_data.dims == 2)
		update_motility_internal<2>(
			data.agents_count, e.mechanics_timestep, data.motility_data.motility_vector.data(), data.velocity.data(),
			data.motility_data.persistence_time.data(), data.motility_data.migration_bias.data(),
			data.motility_data.migration_bias_direction.data(), data.motility_data.restrict_to_2d.data(),
			data.motility_data.is_motile.data(), data.motility_data.migration_speed.data(),
			data.motility_data.direction_update_funcs.data(), data.state_data.agent_type_index.data());
	else if (data.base_data.dims == 3)
		update_motility_internal<3>(
			data.agents_count, e.mechanics_timestep, data.motility_data.motility_vector.data(), data.velocity.data(),
			data.motility_data.persistence_time.data(), data.motility_data.migration_bias.data(),
			data.motility_data.migration_bias_direction.data(), data.motility_data.restrict_to_2d.data(),
			data.motility_data.is_motile.data(), data.motility_data.migration_speed.data(),
			data.motility_data.direction_update_funcs.data(), data.state_data.agent_type_index.data());
}

void position_solver::update_basement_membrane_interactions(environment& e, const cartesian_mesh& mesh)
{
	if (!e.virtual_wall_at_domain_edges) // note: where do we include this
		return;

	auto& data = retrieve_agent_data(*e.agents);

	if (data.base_data.dims == 1)
		update_basement_membrane_interactions_internal<1>(
			data.agents_count, data.velocity.data(), data.base_data.positions.data(), data.radius.data(),
			data.mechanics_data.cell_BM_repulsion_strength.data(), data.state_data.is_movable.data(), mesh);
	else if (data.base_data.dims == 2)
		update_basement_membrane_interactions_internal<2>(
			data.agents_count, data.velocity.data(), data.base_data.positions.data(), data.radius.data(),
			data.mechanics_data.cell_BM_repulsion_strength.data(), data.state_data.is_movable.data(), mesh);
	else if (data.base_data.dims == 3)
		update_basement_membrane_interactions_internal<3>(
			data.agents_count, data.velocity.data(), data.base_data.positions.data(), data.radius.data(),
			data.mechanics_data.cell_BM_repulsion_strength.data(), data.state_data.is_movable.data(), mesh);
}

void position_solver::update_spring_attachments(environment& e)
{
	if (!e.automated_spring_adhesion)
		return;

	auto& data = retrieve_agent_data(*e.agents);
	const index_t dims = data.base_data.dims;

	update_spring_attachments_internal(
		data.agents_count, e.mechanics_timestep, data.agent_types_count, data.mechanics_data.detachment_rate.data(),
		data.mechanics_data.attachment_rate.data(), data.mechanics_data.cell_adhesion_affinities.data(),
		data.mechanics_data.maximum_number_of_attachments.data(), data.state_data.agent_type_index.data(),
		data.state_data.neighbors.data(), data.state_data.springs.data());

	if (dims == 1)
		spring_contract_function<1>(
			data.agents_count, data.agent_types_count, data.velocity.data(), data.state_data.agent_type_index.data(),
			data.mechanics_data.attachment_elastic_constant.data(), data.mechanics_data.cell_adhesion_affinities.data(),
			data.base_data.positions.data(), data.state_data.is_movable.data(), data.state_data.springs.data());
	else if (dims == 2)
		spring_contract_function<2>(
			data.agents_count, data.agent_types_count, data.velocity.data(), data.state_data.agent_type_index.data(),
			data.mechanics_data.attachment_elastic_constant.data(), data.mechanics_data.cell_adhesion_affinities.data(),
			data.base_data.positions.data(), data.state_data.is_movable.data(), data.state_data.springs.data());
	else if (dims == 3)
		spring_contract_function<3>(
			data.agents_count, data.agent_types_count, data.velocity.data(), data.state_data.agent_type_index.data(),
			data.mechanics_data.attachment_elastic_constant.data(), data.mechanics_data.cell_adhesion_affinities.data(),
			data.base_data.positions.data(), data.state_data.is_movable.data(), data.state_data.springs.data());
}

void position_solver::update_positions(environment& e)
{
	auto& data = retrieve_agent_data(*e.agents);
	const index_t dims = data.base_data.dims;

	if (dims == 1)
		update_positions_internal<1>(data.agents_count, e.mechanics_timestep, data.base_data.positions.data(),
									 data.velocity.data(), data.previous_velocity.data(),
									 data.state_data.is_movable.data());
	else if (dims == 2)
		update_positions_internal<2>(data.agents_count, e.mechanics_timestep, data.base_data.positions.data(),
									 data.velocity.data(), data.previous_velocity.data(),
									 data.state_data.is_movable.data());
	else if (dims == 3)
		update_positions_internal<3>(data.agents_count, e.mechanics_timestep, data.base_data.positions.data(),
									 data.velocity.data(), data.previous_velocity.data(),
									 data.state_data.is_movable.data());
}

} // namespace physicore::mechanics::physicell::kernels::openmp_solver
