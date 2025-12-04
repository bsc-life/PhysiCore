#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

#include <common/base_agent_data.h>
#include <common/types.h>

namespace physicore::mechanics::micromechanics {

template <template <typename...> typename ContainerType = std::vector>
struct agent_data_generic_storage
{
public:
	physicore::base_agent_data_generic_storage<ContainerType>& base_data;

	// Agent Classification
	ContainerType<std::uint8_t> agent_types; // Agent type for type-based interactions (0-255)
	ContainerType<index_t> cell_ids;		 // Cell ID this agent belongs to (-1 if standalone)

	// Physics State
	ContainerType<real_t> velocities;		   // n * dims
	ContainerType<real_t> previous_velocities; // n * dims (for Adams-Bashforth)
	ContainerType<real_t> forces;			   // n * dims

	// Geometry & Properties
	ContainerType<real_t> radii;
	ContainerType<std::uint8_t> is_movable;

	// Mechanics Parameters (Per-agent to allow heterogeneity)
	ContainerType<real_t> cell_cell_adhesion_strength;
	ContainerType<real_t> cell_cell_repulsion_strength;
	ContainerType<real_t> relative_maximum_adhesion_distance;

	// Mechanics - Adhesion/Repulsion (Additional)
	ContainerType<real_t> cell_BM_adhesion_strength;
	ContainerType<real_t> cell_BM_repulsion_strength;

	// Mechanics - Attachments
	ContainerType<index_t> maximum_number_of_attachments;
	ContainerType<real_t> attachment_elastic_constant;
	ContainerType<real_t> attachment_rate;
	ContainerType<real_t> detachment_rate;

	// Mechanics - Advanced (Morse / Kelvin-Voigt)
	ContainerType<index_t> cell_residency;
	ContainerType<real_t> intra_scaling_factors;
	ContainerType<real_t> intra_equilibrium_distances;
	ContainerType<real_t> intra_stiffnesses;
	ContainerType<real_t> spring_constants;
	ContainerType<real_t> dissipation_rates;

	// Topology (Kelvin-Voigt)
	// Note: Nested vectors are not ideal for SoA but required for dynamic topology
	ContainerType<std::vector<index_t>> neighbors;
	ContainerType<std::vector<real_t>> rest_lengths;

	// Motility
	ContainerType<std::uint8_t> is_motile;
	ContainerType<real_t> persistence_times;
	ContainerType<real_t> migration_speeds;
	ContainerType<real_t> migration_bias_directions; // n * dims
	ContainerType<real_t> migration_biases;
	ContainerType<real_t> motility_directions; // n * dims (current motility direction)
	ContainerType<std::uint8_t> restrict_to_2d;
	ContainerType<index_t> chemotaxis_index;
	ContainerType<index_t> chemotaxis_direction;

	// Spring attachments (list of attached agent indices per agent)
	ContainerType<std::vector<index_t>> spring_attachments;

	// TODO: Add cell_adhesion_affinities (requires cell_definitions_count)
	// TODO: Add chemotactic_sensitivities (requires substrates_count)



	index_t agents_count = 0;

	explicit agent_data_generic_storage(physicore::base_agent_data_generic_storage<ContainerType>& base_data);

	void add();
	void remove_at(index_t position);
};

template <template <typename...> typename ContainerType>
agent_data_generic_storage<ContainerType>::agent_data_generic_storage(
	physicore::base_agent_data_generic_storage<ContainerType>& base_data)
	: base_data(base_data)
{}

template <template <typename...> typename ContainerType>
void agent_data_generic_storage<ContainerType>::add()
{
	++agents_count;
	index_t dims = base_data.dims;

	// Agent Classification
	agent_types.resize(agents_count, 0);
	cell_ids.resize(agents_count, static_cast<index_t>(-1)); // -1 = standalone agent

	velocities.resize(agents_count * dims, 0.0);
	previous_velocities.resize(agents_count * dims, 0.0);
	forces.resize(agents_count * dims, 0.0);

	radii.resize(agents_count, 0.0);
	is_movable.resize(agents_count, 1);

	cell_cell_adhesion_strength.resize(agents_count, 0.0);
	cell_cell_repulsion_strength.resize(agents_count, 0.0);
	relative_maximum_adhesion_distance.resize(agents_count, 0.0);

	cell_BM_adhesion_strength.resize(agents_count, 0.0);
	cell_BM_repulsion_strength.resize(agents_count, 0.0);

	maximum_number_of_attachments.resize(agents_count, 0);
	attachment_elastic_constant.resize(agents_count, 0.0);
	attachment_rate.resize(agents_count, 0.0);
	detachment_rate.resize(agents_count, 0.0);

	cell_residency.resize(agents_count, 0);
	intra_scaling_factors.resize(agents_count, 0.0);
	intra_equilibrium_distances.resize(agents_count, 0.0);
	intra_stiffnesses.resize(agents_count, 0.0);
	spring_constants.resize(agents_count, 0.0);
	dissipation_rates.resize(agents_count, 0.0);

	neighbors.resize(agents_count);
	rest_lengths.resize(agents_count);

	is_motile.resize(agents_count, 0);
	persistence_times.resize(agents_count, 0.0);
	migration_speeds.resize(agents_count, 0.0);
	migration_bias_directions.resize(agents_count * dims, 0.0);
	migration_biases.resize(agents_count, 0.0);
	motility_directions.resize(agents_count * dims, 0.0);
	restrict_to_2d.resize(agents_count, 0);
	chemotaxis_index.resize(agents_count, 0);
	chemotaxis_direction.resize(agents_count, 0);

	spring_attachments.resize(agents_count);
}

template <template <typename...> typename ContainerType>
void agent_data_generic_storage<ContainerType>::remove_at(index_t position)
{
	assert(position < agents_count);
	if (position >= agents_count)
		return;

	--agents_count;
	index_t dims = base_data.dims;

	if (position < agents_count)
	{
		// Agent Classification
		base_agent_data::move_scalar(&agent_types[position], &agent_types[agents_count]);
		base_agent_data::move_scalar(&cell_ids[position], &cell_ids[agents_count]);

		base_agent_data::move_vector(&velocities[position * dims], &velocities[agents_count * dims], dims);
		base_agent_data::move_vector(&previous_velocities[position * dims], &previous_velocities[agents_count * dims],
									 dims);
		base_agent_data::move_vector(&forces[position * dims], &forces[agents_count * dims], dims);

		base_agent_data::move_scalar(&radii[position], &radii[agents_count]);
		base_agent_data::move_scalar(&is_movable[position], &is_movable[agents_count]);

		base_agent_data::move_scalar(&cell_cell_adhesion_strength[position],
									 &cell_cell_adhesion_strength[agents_count]);
		base_agent_data::move_scalar(&cell_cell_repulsion_strength[position],
									 &cell_cell_repulsion_strength[agents_count]);
		base_agent_data::move_scalar(&relative_maximum_adhesion_distance[position],
									 &relative_maximum_adhesion_distance[agents_count]);

		base_agent_data::move_scalar(&cell_BM_adhesion_strength[position], &cell_BM_adhesion_strength[agents_count]);
		base_agent_data::move_scalar(&cell_BM_repulsion_strength[position], &cell_BM_repulsion_strength[agents_count]);

		base_agent_data::move_scalar(&maximum_number_of_attachments[position],
									 &maximum_number_of_attachments[agents_count]);
		base_agent_data::move_scalar(&attachment_elastic_constant[position],
									 &attachment_elastic_constant[agents_count]);
		base_agent_data::move_scalar(&attachment_rate[position], &attachment_rate[agents_count]);
		base_agent_data::move_scalar(&detachment_rate[position], &detachment_rate[agents_count]);

		base_agent_data::move_scalar(&cell_residency[position], &cell_residency[agents_count]);
		base_agent_data::move_scalar(&intra_scaling_factors[position], &intra_scaling_factors[agents_count]);
		base_agent_data::move_scalar(&intra_equilibrium_distances[position],
									 &intra_equilibrium_distances[agents_count]);
		base_agent_data::move_scalar(&intra_stiffnesses[position], &intra_stiffnesses[agents_count]);
		base_agent_data::move_scalar(&spring_constants[position], &spring_constants[agents_count]);
		base_agent_data::move_scalar(&dissipation_rates[position], &dissipation_rates[agents_count]);

		// Move complex types manually
		neighbors[position] = std::move(neighbors[agents_count]);
		rest_lengths[position] = std::move(rest_lengths[agents_count]);

		base_agent_data::move_scalar(&is_motile[position], &is_motile[agents_count]);
		base_agent_data::move_scalar(&persistence_times[position], &persistence_times[agents_count]);
		base_agent_data::move_scalar(&migration_speeds[position], &migration_speeds[agents_count]);
		base_agent_data::move_vector(&migration_bias_directions[position * dims],
									 &migration_bias_directions[agents_count * dims], dims);
		base_agent_data::move_scalar(&migration_biases[position], &migration_biases[agents_count]);
		base_agent_data::move_vector(&motility_directions[position * dims], &motility_directions[agents_count * dims],
									 dims);
		base_agent_data::move_scalar(&restrict_to_2d[position], &restrict_to_2d[agents_count]);
		base_agent_data::move_scalar(&chemotaxis_index[position], &chemotaxis_index[agents_count]);
		base_agent_data::move_scalar(&chemotaxis_direction[position], &chemotaxis_direction[agents_count]);

		spring_attachments[position] = std::move(spring_attachments[agents_count]);
	}

	// Resize to shrink
	agent_types.resize(agents_count);
	cell_ids.resize(agents_count);
	velocities.resize(agents_count * dims);
	previous_velocities.resize(agents_count * dims);
	forces.resize(agents_count * dims);
	radii.resize(agents_count);
	is_movable.resize(agents_count);
	cell_cell_adhesion_strength.resize(agents_count);
	cell_cell_repulsion_strength.resize(agents_count);
	relative_maximum_adhesion_distance.resize(agents_count);

	cell_BM_adhesion_strength.resize(agents_count);
	cell_BM_repulsion_strength.resize(agents_count);

	maximum_number_of_attachments.resize(agents_count);
	attachment_elastic_constant.resize(agents_count);
	attachment_rate.resize(agents_count);
	detachment_rate.resize(agents_count);

	cell_residency.resize(agents_count);
	intra_scaling_factors.resize(agents_count);
	intra_equilibrium_distances.resize(agents_count);
	intra_stiffnesses.resize(agents_count);
	spring_constants.resize(agents_count);
	dissipation_rates.resize(agents_count);

	neighbors.resize(agents_count);
	rest_lengths.resize(agents_count);

	is_motile.resize(agents_count);
	persistence_times.resize(agents_count);
	migration_speeds.resize(agents_count);
	migration_bias_directions.resize(agents_count * dims);
	migration_biases.resize(agents_count);
	motility_directions.resize(agents_count * dims);
	restrict_to_2d.resize(agents_count);
	chemotaxis_index.resize(agents_count);
	chemotaxis_direction.resize(agents_count);

	spring_attachments.resize(agents_count);
}

using agent_data = agent_data_generic_storage<std::vector>;

} // namespace physicore::mechanics::micromechanics
