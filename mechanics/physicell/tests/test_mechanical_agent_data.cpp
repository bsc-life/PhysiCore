#include <common/base_agent_data.h>
#include <gtest/gtest.h>

#include "physicell/agent_data.h"

namespace physicore::mechanics::physicell::tests {

using namespace physicore;
using namespace physicore::mechanics::physicell;

namespace {
base_agent_data make_base_agent_data(index_t dims, index_t count)
{
	base_agent_data data(dims);
	for (index_t i = 0; i < count; ++i)
	{
		data.add();
	}
	return data;
}

void sync_add(base_agent_data& base, agent_data& data)
{
	base.add();
	data.add();
}

void sync_remove_at(base_agent_data& base, agent_data& data, index_t position)
{
	base.remove_at(position);
	data.remove_at(position);
}

void expect_container_sizes(const agent_data& data, index_t dims, index_t agent_types_count, index_t substrates_count)
{
	const index_t n = data.agents_count;

	// Kinematics
	EXPECT_EQ(data.velocity.size(), n * dims);
	EXPECT_EQ(data.previous_velocity.size(), n * dims);
	EXPECT_EQ(data.radius.size(), n);

	// Mechanics properties
	EXPECT_EQ(data.mechanics_data.cell_cell_adhesion_strength.size(), n);
	EXPECT_EQ(data.mechanics_data.cell_BM_adhesion_strength.size(), n);
	EXPECT_EQ(data.mechanics_data.cell_cell_repulsion_strength.size(), n);
	EXPECT_EQ(data.mechanics_data.cell_BM_repulsion_strength.size(), n);
	EXPECT_EQ(data.mechanics_data.cell_adhesion_affinities.size(), n * agent_types_count);
	EXPECT_EQ(data.mechanics_data.relative_maximum_adhesion_distance.size(), n);
	EXPECT_EQ(data.mechanics_data.maximum_number_of_attachments.size(), n);
	EXPECT_EQ(data.mechanics_data.attachment_elastic_constant.size(), n);
	EXPECT_EQ(data.mechanics_data.attachment_rate.size(), n);
	EXPECT_EQ(data.mechanics_data.detachment_rate.size(), n);

	// Motility properties
	EXPECT_EQ(data.motility_data.is_motile.size(), n);
	EXPECT_EQ(data.motility_data.persistence_time.size(), n);
	EXPECT_EQ(data.motility_data.migration_speed.size(), n);
	EXPECT_EQ(data.motility_data.migration_bias_direction.size(), n * dims);
	EXPECT_EQ(data.motility_data.migration_bias.size(), n);
	EXPECT_EQ(data.motility_data.motility_vector.size(), n * dims);
	EXPECT_EQ(data.motility_data.restrict_to_2d.size(), n);
	EXPECT_EQ(data.motility_data.chemotaxis_index.size(), n);
	EXPECT_EQ(data.motility_data.chemotaxis_direction.size(), n);
	EXPECT_EQ(data.motility_data.chemotactic_sensitivities.size(), n * substrates_count);

	// State properties
	EXPECT_EQ(data.state_data.neighbors.size(), n);
	EXPECT_EQ(data.state_data.springs.size(), n);
	EXPECT_EQ(data.state_data.attached_cells.size(), n);
	EXPECT_EQ(data.state_data.orientation.size(), n * dims);
	EXPECT_EQ(data.state_data.simple_pressure.size(), n);
	EXPECT_EQ(data.state_data.agent_type_index.size(), n);
	EXPECT_EQ(data.state_data.is_movable.size(), n);
}
} // namespace

TEST(MechanicalAgentDataTest, AddInitializesVectorsCorrectly)
{
	const index_t dims = 3;
	base_agent_data base = make_base_agent_data(dims, 0);
	const index_t agent_types_count = 4;
	const index_t substrates_count = 2;
	agent_data data(base, agent_types_count, substrates_count);

	EXPECT_EQ(data.agents_count, 0);
	EXPECT_EQ(data.base_data.dims, dims);
	EXPECT_EQ(data.agent_types_count, agent_types_count);
	EXPECT_EQ(data.substrates_count, substrates_count);
	expect_container_sizes(data, dims, agent_types_count, substrates_count);

	sync_add(base, data);

	EXPECT_EQ(data.agents_count, 1);
	EXPECT_EQ(data.agents_count, base.agents_count);
	expect_container_sizes(data, dims, agent_types_count, substrates_count);

	sync_add(base, data);

	EXPECT_EQ(data.agents_count, 2);
	EXPECT_EQ(data.agents_count, base.agents_count);
	expect_container_sizes(data, dims, agent_types_count, substrates_count);
}

TEST(MechanicalAgentDataTest, RemoveShrinksVectorsCorrectly)
{
	const index_t dims = 3;
	base_agent_data base = make_base_agent_data(dims, 0);
	const index_t agent_types_count = 3;
	const index_t substrates_count = 2;
	agent_data data(base, agent_types_count, substrates_count);

	sync_add(base, data); // agent 0
	sync_add(base, data); // agent 1
	sync_add(base, data); // agent 2

	EXPECT_EQ(data.agents_count, 3);
	EXPECT_EQ(data.agents_count, base.agents_count);
	expect_container_sizes(data, dims, agent_types_count, substrates_count);

	// Remove agent at position 1
	sync_remove_at(base, data, 1);

	EXPECT_EQ(data.agents_count, 2);
	EXPECT_EQ(data.agents_count, base.agents_count);
	expect_container_sizes(data, dims, agent_types_count, substrates_count);

	// Remove agent at position 0
	sync_remove_at(base, data, 0);

	EXPECT_EQ(data.agents_count, 1);
	EXPECT_EQ(data.agents_count, base.agents_count);
	expect_container_sizes(data, dims, agent_types_count, substrates_count);
}

TEST(MechanicalAgentDataTest, RemoveMovesLastAgentDataToRemovedSlot)
{
	const index_t dims = 3;
	base_agent_data base = make_base_agent_data(dims, 0);
	const index_t agent_types_count = 2;
	const index_t substrates_count = 2;
	agent_data data(base, agent_types_count, substrates_count);

	sync_add(base, data); // agent 0
	sync_add(base, data); // agent 1
	sync_add(base, data); // agent 2

	// Stamp distinct values in multiple containers.
	data.radius[0] = 10.0;
	data.radius[1] = 20.0;
	data.radius[2] = 30.0;

	data.velocity[0 * dims + 0] = 0.1;
	data.velocity[0 * dims + 1] = 0.2;
	data.velocity[0 * dims + 2] = 0.3;
	data.velocity[1 * dims + 0] = 1.1;
	data.velocity[1 * dims + 1] = 1.2;
	data.velocity[1 * dims + 2] = 1.3;
	data.velocity[2 * dims + 0] = 2.1;
	data.velocity[2 * dims + 1] = 2.2;
	data.velocity[2 * dims + 2] = 2.3;

	data.state_data.orientation[0 * dims + 0] = 10.0;
	data.state_data.orientation[0 * dims + 1] = 20.0;
	data.state_data.orientation[0 * dims + 2] = 30.0;
	data.state_data.orientation[1 * dims + 0] = 11.0;
	data.state_data.orientation[1 * dims + 1] = 21.0;
	data.state_data.orientation[1 * dims + 2] = 31.0;
	data.state_data.orientation[2 * dims + 0] = 12.0;
	data.state_data.orientation[2 * dims + 1] = 22.0;
	data.state_data.orientation[2 * dims + 2] = 32.0;

	data.mechanics_data.cell_adhesion_affinities[0 * agent_types_count + 0] = 0.01;
	data.mechanics_data.cell_adhesion_affinities[0 * agent_types_count + 1] = 0.02;
	data.mechanics_data.cell_adhesion_affinities[1 * agent_types_count + 0] = 0.11;
	data.mechanics_data.cell_adhesion_affinities[1 * agent_types_count + 1] = 0.12;
	data.mechanics_data.cell_adhesion_affinities[2 * agent_types_count + 0] = 0.21;
	data.mechanics_data.cell_adhesion_affinities[2 * agent_types_count + 1] = 0.22;

	data.motility_data.chemotactic_sensitivities[0 * substrates_count + 0] = 1.0;
	data.motility_data.chemotactic_sensitivities[0 * substrates_count + 1] = 2.0;
	data.motility_data.chemotactic_sensitivities[1 * substrates_count + 0] = 3.0;
	data.motility_data.chemotactic_sensitivities[1 * substrates_count + 1] = 4.0;
	data.motility_data.chemotactic_sensitivities[2 * substrates_count + 0] = 5.0;
	data.motility_data.chemotactic_sensitivities[2 * substrates_count + 1] = 6.0;

	data.state_data.neighbors[0] = { 1 };
	data.state_data.neighbors[1] = { 0, 2 };
	data.state_data.neighbors[2] = { 42 };

	// Remove the middle agent, expect last data to move into slot 1.
	sync_remove_at(base, data, 1);

	EXPECT_EQ(data.agents_count, 2);
	EXPECT_EQ(data.agents_count, base.agents_count);

	// Slot 0 unchanged.
	EXPECT_DOUBLE_EQ(data.radius[0], 10.0);
	EXPECT_DOUBLE_EQ(data.velocity[0 * dims + 0], 0.1);
	EXPECT_DOUBLE_EQ(data.velocity[0 * dims + 1], 0.2);
	EXPECT_DOUBLE_EQ(data.velocity[0 * dims + 2], 0.3);
	EXPECT_DOUBLE_EQ(data.state_data.orientation[0 * dims + 0], 10.0);
	EXPECT_DOUBLE_EQ(data.state_data.orientation[0 * dims + 1], 20.0);
	EXPECT_DOUBLE_EQ(data.state_data.orientation[0 * dims + 2], 30.0);
	ASSERT_EQ(data.state_data.neighbors[0].size(), 1);
	EXPECT_EQ(data.state_data.neighbors[0][0], 1);

	// Slot 1 now contains what was previously slot 2.
	EXPECT_DOUBLE_EQ(data.radius[1], 30.0);
	EXPECT_DOUBLE_EQ(data.velocity[1 * dims + 0], 2.1);
	EXPECT_DOUBLE_EQ(data.velocity[1 * dims + 1], 2.2);
	EXPECT_DOUBLE_EQ(data.velocity[1 * dims + 2], 2.3);
	EXPECT_DOUBLE_EQ(data.state_data.orientation[1 * dims + 0], 12.0);
	EXPECT_DOUBLE_EQ(data.state_data.orientation[1 * dims + 1], 22.0);
	EXPECT_DOUBLE_EQ(data.state_data.orientation[1 * dims + 2], 32.0);

	EXPECT_DOUBLE_EQ(data.mechanics_data.cell_adhesion_affinities[1 * agent_types_count + 0], 0.21);
	EXPECT_DOUBLE_EQ(data.mechanics_data.cell_adhesion_affinities[1 * agent_types_count + 1], 0.22);
	EXPECT_DOUBLE_EQ(data.motility_data.chemotactic_sensitivities[1 * substrates_count + 0], 5.0);
	EXPECT_DOUBLE_EQ(data.motility_data.chemotactic_sensitivities[1 * substrates_count + 1], 6.0);

	ASSERT_EQ(data.state_data.neighbors[1].size(), 1);
	EXPECT_EQ(data.state_data.neighbors[1][0], 42);
}

} // namespace physicore::mechanics::physicell::tests
