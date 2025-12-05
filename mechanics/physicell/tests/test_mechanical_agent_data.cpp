#include <memory>
#include <vector>

#include <common/base_agent_data.h>
#include <common/types.h>
#include <gtest/gtest.h>

#include "physicell/agent_data.h"

namespace physicore::mechanics::physicell::tests {

using namespace physicore;
using namespace physicore::mechanics::physicell;

class MechanicalAgentDataTest : public ::testing::Test
{
protected:
	std::unique_ptr<base_agent_data_generic_storage<std::vector>> base_data;

	void SetUp() override { base_data = std::make_unique<base_agent_data_generic_storage<std::vector>>(3); }
};

// Constructor tests
TEST_F(MechanicalAgentDataTest, ConstructorDefaultValues)
{
	mechanical_agent_data agent_data(*base_data);

	EXPECT_EQ(agent_data.agents_count, 0);
	EXPECT_EQ(agent_data.dims, 3);
	EXPECT_EQ(agent_data.agent_types_count, 0);
	EXPECT_EQ(agent_data.substrates_count, 0);
}

TEST_F(MechanicalAgentDataTest, ConstructorCustomDimensions)
{
	// dims is taken from base_data.dims if > 0, otherwise from parameter
	mechanical_agent_data agent_data(*base_data, 2, 5, 3);

	EXPECT_EQ(agent_data.dims, 3); // Takes base_data.dims (3)
	EXPECT_EQ(agent_data.agent_types_count, 5);
	EXPECT_EQ(agent_data.substrates_count, 3);
	EXPECT_EQ(agent_data.agents_count, 0);
}

TEST_F(MechanicalAgentDataTest, ConstructorWithDifferentParams)
{
	base_data = std::make_unique<base_agent_data_generic_storage<std::vector>>(3);
	mechanical_agent_data agent_data(*base_data, 3, 10, 5);

	// dims comes from base_data
	EXPECT_EQ(agent_data.dims, 3);
	EXPECT_EQ(agent_data.agent_types_count, 10);
	EXPECT_EQ(agent_data.substrates_count, 5);
}

// Add single agent tests
TEST_F(MechanicalAgentDataTest, AddSingleAgent)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 1);

	base_data->add();
	agent_data.add();

	EXPECT_EQ(agent_data.agents_count, 1);
	EXPECT_EQ(agent_data.velocities.size(), 3); // dims
	EXPECT_EQ(agent_data.previous_velocities.size(), 3);
	EXPECT_EQ(agent_data.radius.size(), 1);
}

TEST_F(MechanicalAgentDataTest, AddMultipleAgents)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 1);

	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();

	EXPECT_EQ(agent_data.agents_count, 3);
	EXPECT_EQ(agent_data.velocities.size(), 9); // 3 agents * 3 dims
	EXPECT_EQ(agent_data.radius.size(), 3);
	EXPECT_EQ(agent_data.is_motile.size(), 3);
}

// Adhesion container tests
TEST_F(MechanicalAgentDataTest, AdhesionContainersAfterAdd)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 1);

	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();

	EXPECT_EQ(agent_data.cell_cell_adhesion_strength.size(), 2);
	EXPECT_EQ(agent_data.cell_BM_adhesion_strength.size(), 2);
	EXPECT_EQ(agent_data.cell_cell_repulsion_strength.size(), 2);
	EXPECT_EQ(agent_data.cell_BM_repulsion_strength.size(), 2);
	EXPECT_EQ(agent_data.relative_maximum_adhesion_distance.size(), 2);
}

// Motility container tests
TEST_F(MechanicalAgentDataTest, MotilityContainersAfterAdd)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 1);

	base_data->add();
	agent_data.add();

	EXPECT_EQ(agent_data.is_motile.size(), 1);
	EXPECT_EQ(agent_data.persistence_time.size(), 1);
	EXPECT_EQ(agent_data.migration_speed.size(), 1);
	EXPECT_EQ(agent_data.migration_bias.size(), 1);
	EXPECT_EQ(agent_data.restrict_to_2d.size(), 1);
}

// Chemotaxis container tests
TEST_F(MechanicalAgentDataTest, ChemotaxisContainersAfterAdd)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 2);

	base_data->add();
	agent_data.add();

	EXPECT_EQ(agent_data.chemotaxis_index.size(), 1);
	EXPECT_EQ(agent_data.chemotaxis_direction.size(), 1);
	EXPECT_EQ(agent_data.chemotactic_sensitivities.size(), 2); // substrates_count
}

// Affinities container test
TEST_F(MechanicalAgentDataTest, AffinitiiesContainerAfterAdd)
{
	mechanical_agent_data agent_data(*base_data, 3, 3, 1);

	base_data->add();
	agent_data.add();

	EXPECT_EQ(agent_data.cell_adhesion_affinities.size(), 3); // agent_types_count
}

// State container tests
TEST_F(MechanicalAgentDataTest, StateContainersAfterAdd)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 1);

	base_data->add();
	agent_data.add();

	EXPECT_EQ(agent_data.neighbors.size(), 1);
	EXPECT_EQ(agent_data.springs.size(), 1);
	EXPECT_EQ(agent_data.attached_cells.size(), 1);
	EXPECT_EQ(agent_data.orientation.size(), 3); // dims
	EXPECT_EQ(agent_data.simple_pressure.size(), 1);
	EXPECT_EQ(agent_data.cell_definition_indices.size(), 1);
	EXPECT_EQ(agent_data.is_movable.size(), 1);
}

// Remove tests
TEST_F(MechanicalAgentDataTest, RemoveAtFirstAgent)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 1);

	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();
	EXPECT_EQ(agent_data.agents_count, 3);

	base_data->remove_at(0);
	agent_data.remove_at(0);

	EXPECT_EQ(agent_data.agents_count, 2);
	EXPECT_EQ(agent_data.velocities.size(), 6); // 2 agents * 3 dims
	EXPECT_EQ(agent_data.radius.size(), 2);
}

TEST_F(MechanicalAgentDataTest, RemoveAtLastAgent)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 1);

	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();

	base_data->remove_at(2);
	agent_data.remove_at(2);

	EXPECT_EQ(agent_data.agents_count, 2);
	EXPECT_EQ(agent_data.velocities.size(), 6);
}

TEST_F(MechanicalAgentDataTest, RemoveAtMiddleAgent)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 1);

	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();

	base_data->remove_at(1);
	agent_data.remove_at(1);

	EXPECT_EQ(agent_data.agents_count, 2);
	EXPECT_EQ(agent_data.velocities.size(), 6);
	EXPECT_EQ(agent_data.radius.size(), 2);
}

TEST_F(MechanicalAgentDataTest, RemoveSingleAgent)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 1);

	base_data->add();
	agent_data.add();
	EXPECT_EQ(agent_data.agents_count, 1);

	base_data->remove_at(0);
	agent_data.remove_at(0);

	EXPECT_EQ(agent_data.agents_count, 0);
	EXPECT_EQ(agent_data.velocities.size(), 0);
}

// Multiple add/remove cycles
TEST_F(MechanicalAgentDataTest, MultipleAddRemoveCycles)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 1);

	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();
	EXPECT_EQ(agent_data.agents_count, 2);

	base_data->remove_at(0);
	agent_data.remove_at(0);
	EXPECT_EQ(agent_data.agents_count, 1);

	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();
	EXPECT_EQ(agent_data.agents_count, 3);

	base_data->remove_at(1);
	agent_data.remove_at(1);
	EXPECT_EQ(agent_data.agents_count, 2);
}

// Container consistency tests
TEST_F(MechanicalAgentDataTest, ContainerConsistencyAfterAddRemove)
{
	mechanical_agent_data agent_data(*base_data, 3, 4, 2);

	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();

	// All scalar containers should have size = agents_count
	EXPECT_EQ(agent_data.velocities.size(), agent_data.agents_count * 3);
	EXPECT_EQ(agent_data.radius.size(), agent_data.agents_count);
	EXPECT_EQ(agent_data.is_motile.size(), agent_data.agents_count);
	EXPECT_EQ(agent_data.is_movable.size(), agent_data.agents_count);

	base_data->remove_at(1);
	agent_data.remove_at(1);

	// Consistency after removal
	EXPECT_EQ(agent_data.velocities.size(), agent_data.agents_count * 3);
	EXPECT_EQ(agent_data.radius.size(), agent_data.agents_count);
	EXPECT_EQ(agent_data.is_motile.size(), agent_data.agents_count);
	EXPECT_EQ(agent_data.is_movable.size(), agent_data.agents_count);
}

// Attachment container tests
TEST_F(MechanicalAgentDataTest, AttachmentContainersAfterAdd)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 1);

	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();

	EXPECT_EQ(agent_data.maximum_number_of_attachments.size(), 2);
	EXPECT_EQ(agent_data.attachment_elastic_constant.size(), 2);
	EXPECT_EQ(agent_data.attachment_rate.size(), 2);
	EXPECT_EQ(agent_data.detachment_rate.size(), 2);
}

// Dimension-dependent container tests
TEST_F(MechanicalAgentDataTest, DimensionDependentContainersAfterAdd2D)
{
	base_data = std::make_unique<base_agent_data_generic_storage<std::vector>>(2);
	mechanical_agent_data agent_data(*base_data, 2, 1, 1);

	base_data->add();
	agent_data.add();

	// Containers with dims scaling
	EXPECT_EQ(agent_data.velocities.size(), 2); // dims
	EXPECT_EQ(agent_data.previous_velocities.size(), 2);
	EXPECT_EQ(agent_data.migration_bias_direction.size(), 2);
	EXPECT_EQ(agent_data.motility_vector.size(), 2);
	EXPECT_EQ(agent_data.orientation.size(), 2);
}

TEST_F(MechanicalAgentDataTest, DimensionDependentContainersAfterAdd3D)
{
	base_data = std::make_unique<base_agent_data_generic_storage<std::vector>>(3);
	mechanical_agent_data agent_data(*base_data, 3, 1, 1);

	base_data->add();
	agent_data.add();

	EXPECT_EQ(agent_data.velocities.size(), 3);
	EXPECT_EQ(agent_data.previous_velocities.size(), 3);
	EXPECT_EQ(agent_data.migration_bias_direction.size(), 3);
	EXPECT_EQ(agent_data.motility_vector.size(), 3);
	EXPECT_EQ(agent_data.orientation.size(), 3);
}

// Type-dependent container tests
TEST_F(MechanicalAgentDataTest, TypeDependentContainersAfterAdd)
{
	mechanical_agent_data agent_data(*base_data, 3, 5, 1);

	base_data->add();
	agent_data.add();

	EXPECT_EQ(agent_data.cell_adhesion_affinities.size(), 5);
}

TEST_F(MechanicalAgentDataTest, TypeDependentContainersMultipleAgents)
{
	mechanical_agent_data agent_data(*base_data, 3, 4, 1);

	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();

	EXPECT_EQ(agent_data.cell_adhesion_affinities.size(), 8); // 2 agents * 4 types
}

// Substrate-dependent container tests
TEST_F(MechanicalAgentDataTest, SubstrateDependentContainersAfterAdd)
{
	mechanical_agent_data agent_data(*base_data, 3, 1, 4);

	base_data->add();
	agent_data.add();

	EXPECT_EQ(agent_data.chemotactic_sensitivities.size(), 4);
}

TEST_F(MechanicalAgentDataTest, SubstrateDependentContainersMultipleAgents)
{
	mechanical_agent_data agent_data(*base_data, 3, 1, 3);

	base_data->add();
	agent_data.add();
	base_data->add();
	agent_data.add();

	EXPECT_EQ(agent_data.chemotactic_sensitivities.size(), 6); // 2 agents * 3 substrates
}

// Large-scale tests
TEST_F(MechanicalAgentDataTest, LargeScaleAdd)
{
	mechanical_agent_data agent_data(*base_data, 3, 10, 5);

	for (int i = 0; i < 100; ++i)
	{
		base_data->add();
		agent_data.add();
	}
	EXPECT_EQ(agent_data.agents_count, 100);
	EXPECT_EQ(agent_data.velocities.size(), 300); // 100 * 3
	EXPECT_EQ(agent_data.radius.size(), 100);
}

TEST_F(MechanicalAgentDataTest, LargeScaleRemove)
{
	mechanical_agent_data agent_data(*base_data, 3, 10, 5);

	for (int i = 0; i < 100; ++i)
	{
		base_data->add();
		agent_data.add();
	}
	EXPECT_EQ(agent_data.agents_count, 100);

	for (int i = 0; i < 50; ++i)
	{
		base_data->remove_at(0);
		agent_data.remove_at(0);
	}
	EXPECT_EQ(agent_data.agents_count, 50);

	// Container consistency
	EXPECT_EQ(agent_data.velocities.size(), 150); // 50 * 3
	EXPECT_EQ(agent_data.radius.size(), 50);
	EXPECT_EQ(agent_data.cell_definition_indices.size(), 50);
}

// Edge case tests
TEST_F(MechanicalAgentDataTest, AllContainerTypesPresent)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 2);

	base_data->add();
	agent_data.add();

	// Kinematics
	EXPECT_EQ(agent_data.velocities.size(), 3);
	EXPECT_EQ(agent_data.previous_velocities.size(), 3);
	EXPECT_EQ(agent_data.radius.size(), 1);

	// Mechanics
	EXPECT_EQ(agent_data.cell_cell_adhesion_strength.size(), 1);
	EXPECT_EQ(agent_data.cell_adhesion_affinities.size(), 2);

	// Motility
	EXPECT_EQ(agent_data.is_motile.size(), 1);
	EXPECT_EQ(agent_data.migration_bias_direction.size(), 3);

	// State
	EXPECT_EQ(agent_data.neighbors.size(), 1);
	EXPECT_EQ(agent_data.orientation.size(), 3);
}

TEST_F(MechanicalAgentDataTest, EmptyToPopulatedTransition)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 1);

	EXPECT_EQ(agent_data.agents_count, 0);
	EXPECT_EQ(agent_data.velocities.size(), 0);
	EXPECT_EQ(agent_data.radius.size(), 0);

	base_data->add();
	agent_data.add();

	EXPECT_EQ(agent_data.agents_count, 1);
	EXPECT_EQ(agent_data.velocities.size(), 3);
	EXPECT_EQ(agent_data.radius.size(), 1);
}

TEST_F(MechanicalAgentDataTest, PopulatedToEmptyTransition)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 1);

	base_data->add();
	agent_data.add();
	EXPECT_EQ(agent_data.agents_count, 1);

	base_data->remove_at(0);
	agent_data.remove_at(0);

	EXPECT_EQ(agent_data.agents_count, 0);
	EXPECT_EQ(agent_data.velocities.size(), 0);
	EXPECT_EQ(agent_data.radius.size(), 0);
}

// Data integrity test
TEST_F(MechanicalAgentDataTest, DataIntegrityAfterOperations)
{
	mechanical_agent_data agent_data(*base_data, 3, 2, 1);

	// Add agents and store initial values
	base_data->add();
	agent_data.add();
	agent_data.radius[0] = 5.5F;
	agent_data.is_motile[0] = 1;

	base_data->add();
	agent_data.add();
	agent_data.radius[1] = 7.2F;
	agent_data.is_motile[1] = 0;

	base_data->add();
	agent_data.add();
	agent_data.radius[2] = 6.3F;
	agent_data.is_motile[2] = 1;

	// Remove middle agent
	base_data->remove_at(1);
	agent_data.remove_at(1);

	// Check that last agent moved to middle position
	EXPECT_EQ(agent_data.agents_count, 2);
	EXPECT_EQ(agent_data.radius[1], 6.3F);
	EXPECT_EQ(agent_data.is_motile[1], 1);

	// Check that first agent is unchanged
	EXPECT_EQ(agent_data.radius[0], 5.5F);
	EXPECT_EQ(agent_data.is_motile[0], 1);
}

} // namespace physicore::mechanics::physicell::tests
