#include <array>

#include <common/base_agent_data.h>
#include <gtest/gtest.h>

#include "physicell/mechanical_agent.h"

namespace physicore::mechanics::physicell::tests {

using namespace physicore;
using namespace physicore::mechanics::physicell;

class MechanicalAgentTest : public ::testing::Test
{
protected:
	static constexpr index_t dims = 3;
	static constexpr index_t agent_types_count = 4;
	static constexpr index_t substrates_count = 2;

	base_agent_data base_data { dims };
	agent_data data { base_data, agent_types_count, substrates_count };

	template <class Getter>
	void set_and_expect_vec3(Getter&& getter, const std::array<real_t, 3>& values)
	{
		auto vec = getter();
		ASSERT_EQ(vec.size(), static_cast<std::size_t>(dims));

		for (std::size_t i = 0; i < values.size(); ++i)
			vec[i] = values[i];

		for (std::size_t i = 0; i < values.size(); ++i)
			EXPECT_DOUBLE_EQ(getter()[i], values[i]);
	}

	void SetUp() override
	{
		// Add first agent
		base_data.add();
		data.add();
	}
};

TEST_F(MechanicalAgentTest, Position)
{
	mechanical_agent agent(0, data);
	set_and_expect_vec3([&]() { return agent.position(); }, { 1.0, 2.0, 3.0 });
}

TEST_F(MechanicalAgentTest, Velocity)
{
	mechanical_agent agent(0, data);
	set_and_expect_vec3([&]() { return agent.velocity(); }, { 0.1, 0.2, 0.3 });
}

TEST_F(MechanicalAgentTest, PreviousVelocity)
{
	mechanical_agent agent(0, data);
	set_and_expect_vec3([&]() { return agent.previous_velocity(); }, { 1.1, 1.2, 1.3 });
}

TEST_F(MechanicalAgentTest, Radius)
{
	mechanical_agent agent(0, data);
	agent.radius() = 7.5;
	EXPECT_DOUBLE_EQ(agent.radius(), 7.5);
}

TEST_F(MechanicalAgentTest, AdhesionStrengths)
{
	mechanical_agent agent(0, data);

	agent.cell_cell_adhesion_strength() = 1.5;
	agent.cell_BM_adhesion_strength() = 2.5;

	EXPECT_DOUBLE_EQ(agent.cell_cell_adhesion_strength(), 1.5);
	EXPECT_DOUBLE_EQ(agent.cell_BM_adhesion_strength(), 2.5);
}

TEST_F(MechanicalAgentTest, RepulsionStrengths)
{
	mechanical_agent agent(0, data);

	agent.cell_cell_repulsion_strength() = 3.5;
	agent.cell_BM_repulsion_strength() = 4.5;

	EXPECT_DOUBLE_EQ(agent.cell_cell_repulsion_strength(), 3.5);
	EXPECT_DOUBLE_EQ(agent.cell_BM_repulsion_strength(), 4.5);
}

TEST_F(MechanicalAgentTest, CellAdhesionAffinities)
{
	mechanical_agent agent(0, data);
	auto affinities = agent.cell_adhesion_affinities();
	ASSERT_EQ(affinities.size(), agent_types_count);

	affinities[0] = 0.01;
	affinities[1] = 0.02;
	affinities[2] = 0.03;
	affinities[3] = 0.04;

	EXPECT_DOUBLE_EQ(agent.cell_adhesion_affinities()[0], 0.01);
	EXPECT_DOUBLE_EQ(agent.cell_adhesion_affinities()[1], 0.02);
	EXPECT_DOUBLE_EQ(agent.cell_adhesion_affinities()[2], 0.03);
	EXPECT_DOUBLE_EQ(agent.cell_adhesion_affinities()[3], 0.04);
}

TEST_F(MechanicalAgentTest, AttachmentParameters)
{
	mechanical_agent agent(0, data);

	agent.relative_maximum_adhesion_distance() = 1.2;
	agent.maximum_number_of_attachments() = 13;
	agent.attachment_elastic_constant() = 0.4;
	agent.attachment_rate() = 0.05;
	agent.detachment_rate() = 0.06;

	EXPECT_DOUBLE_EQ(agent.relative_maximum_adhesion_distance(), 1.2);
	EXPECT_EQ(agent.maximum_number_of_attachments(), static_cast<index_t>(13));
	EXPECT_DOUBLE_EQ(agent.attachment_elastic_constant(), 0.4);
	EXPECT_DOUBLE_EQ(agent.attachment_rate(), 0.05);
	EXPECT_DOUBLE_EQ(agent.detachment_rate(), 0.06);
}

TEST_F(MechanicalAgentTest, MotilityParameters)
{
	mechanical_agent agent(0, data);

	agent.is_motile() = static_cast<std::uint8_t>(1);
	agent.persistence_time() = 2.0;
	agent.migration_speed() = 3.0;
	agent.migration_bias() = 0.7;
	agent.restrict_to_2d() = static_cast<std::uint8_t>(1);

	EXPECT_EQ(agent.is_motile(), static_cast<std::uint8_t>(1));
	EXPECT_DOUBLE_EQ(agent.persistence_time(), 2.0);
	EXPECT_DOUBLE_EQ(agent.migration_speed(), 3.0);
	EXPECT_DOUBLE_EQ(agent.migration_bias(), 0.7);
	EXPECT_EQ(agent.restrict_to_2d(), static_cast<std::uint8_t>(1));
}

TEST_F(MechanicalAgentTest, MigrationBiasDirection)
{
	mechanical_agent agent(0, data);
	set_and_expect_vec3([&]() { return agent.migration_bias_direction(); }, { 0.3, 0.4, 0.5 });
}

TEST_F(MechanicalAgentTest, MotilityVector)
{
	mechanical_agent agent(0, data);
	set_and_expect_vec3([&]() { return agent.motility_vector(); }, { 4.1, 4.2, 4.3 });
}

TEST_F(MechanicalAgentTest, ChemotaxisParameters)
{
	mechanical_agent agent(0, data);

	agent.chemotaxis_index() = 1;
	agent.chemotaxis_direction() = 2;

	auto sensitivities = agent.chemotactic_sensitivities();
	ASSERT_EQ(sensitivities.size(), substrates_count);

	sensitivities[0] = 0.9;
	sensitivities[1] = 1.1;

	EXPECT_EQ(agent.chemotaxis_index(), static_cast<index_t>(1));
	EXPECT_EQ(agent.chemotaxis_direction(), static_cast<index_t>(2));
	EXPECT_DOUBLE_EQ(agent.chemotactic_sensitivities()[0], 0.9);
	EXPECT_DOUBLE_EQ(agent.chemotactic_sensitivities()[1], 1.1);
}

TEST_F(MechanicalAgentTest, NeighborsSpringsAndAttachedCells)
{
	mechanical_agent agent(0, data);
	EXPECT_EQ(agent.neighbors().size(), 0U);
	EXPECT_EQ(agent.springs().size(), 0U);
	EXPECT_EQ(agent.attached_cells().size(), 0U);

	data.state_data.neighbors[0] = { 10, 20, 30 };
	data.state_data.springs[0] = { 1, 2 };
	data.state_data.attached_cells[0] = { 99 };

	auto neighbors = agent.neighbors();
	auto springs = agent.springs();
	auto attached = agent.attached_cells();

	ASSERT_EQ(neighbors.size(), 3U);
	EXPECT_EQ(neighbors[0], static_cast<index_t>(10));
	EXPECT_EQ(neighbors[1], static_cast<index_t>(20));
	EXPECT_EQ(neighbors[2], static_cast<index_t>(30));

	neighbors[1] = 42;
	EXPECT_EQ(data.state_data.neighbors[0][1], static_cast<index_t>(42));

	ASSERT_EQ(springs.size(), 2U);
	EXPECT_EQ(springs[0], static_cast<index_t>(1));
	EXPECT_EQ(springs[1], static_cast<index_t>(2));

	springs[0] = 123;
	EXPECT_EQ(data.state_data.springs[0][0], static_cast<index_t>(123));

	ASSERT_EQ(attached.size(), 1U);
	EXPECT_EQ(attached[0], static_cast<index_t>(99));

	attached[0] = 77;
	EXPECT_EQ(data.state_data.attached_cells[0][0], static_cast<index_t>(77));
}

TEST_F(MechanicalAgentTest, Orientation)
{
	mechanical_agent agent(0, data);
	set_and_expect_vec3([&]() { return agent.orientation(); }, { 10.0, 20.0, 30.0 });
}

TEST_F(MechanicalAgentTest, StateScalars)
{
	mechanical_agent agent(0, data);

	agent.simple_pressure() = 0.8;
	agent.agent_type_index() = 3;
	agent.is_movable() = static_cast<std::uint8_t>(1);

	EXPECT_DOUBLE_EQ(agent.simple_pressure(), 0.8);
	EXPECT_EQ(agent.agent_type_index(), static_cast<index_t>(3));
	EXPECT_EQ(agent.is_movable(), static_cast<std::uint8_t>(1));
}

TEST_F(MechanicalAgentTest, MultipleAgentsAccessDifferentMemoryLocations)
{
	// Add second agent
	base_data.add();
	data.add();

	mechanical_agent agent0(0, data);
	mechanical_agent agent1(1, data);

	agent0.radius() = 5.0;
	agent1.radius() = 6.0;

	agent0.velocity()[0] = 0.1;
	agent0.velocity()[1] = 0.2;
	agent0.velocity()[2] = 0.3;

	agent1.velocity()[0] = 1.1;
	agent1.velocity()[1] = 1.2;
	agent1.velocity()[2] = 1.3;

	agent0.cell_adhesion_affinities()[0] = 0.01;
	agent0.cell_adhesion_affinities()[1] = 0.02;
	agent0.cell_adhesion_affinities()[2] = 0.03;
	agent0.cell_adhesion_affinities()[3] = 0.04;

	agent1.cell_adhesion_affinities()[0] = 0.11;
	agent1.cell_adhesion_affinities()[1] = 0.12;
	agent1.cell_adhesion_affinities()[2] = 0.13;
	agent1.cell_adhesion_affinities()[3] = 0.14;

	agent0.chemotactic_sensitivities()[0] = 2.0;
	agent0.chemotactic_sensitivities()[1] = 3.0;
	agent1.chemotactic_sensitivities()[0] = 4.0;
	agent1.chemotactic_sensitivities()[1] = 5.0;

	data.state_data.neighbors[0] = { 1 };
	data.state_data.neighbors[1] = { 0, 2 };

	EXPECT_DOUBLE_EQ(agent0.radius(), 5.0);
	EXPECT_DOUBLE_EQ(agent1.radius(), 6.0);

	EXPECT_DOUBLE_EQ(agent0.velocity()[0], 0.1);
	EXPECT_DOUBLE_EQ(agent0.velocity()[1], 0.2);
	EXPECT_DOUBLE_EQ(agent0.velocity()[2], 0.3);

	EXPECT_DOUBLE_EQ(agent1.velocity()[0], 1.1);
	EXPECT_DOUBLE_EQ(agent1.velocity()[1], 1.2);
	EXPECT_DOUBLE_EQ(agent1.velocity()[2], 1.3);

	EXPECT_DOUBLE_EQ(agent0.cell_adhesion_affinities()[0], 0.01);
	EXPECT_DOUBLE_EQ(agent0.cell_adhesion_affinities()[1], 0.02);
	EXPECT_DOUBLE_EQ(agent0.cell_adhesion_affinities()[2], 0.03);
	EXPECT_DOUBLE_EQ(agent0.cell_adhesion_affinities()[3], 0.04);

	EXPECT_DOUBLE_EQ(agent1.cell_adhesion_affinities()[0], 0.11);
	EXPECT_DOUBLE_EQ(agent1.cell_adhesion_affinities()[1], 0.12);
	EXPECT_DOUBLE_EQ(agent1.cell_adhesion_affinities()[2], 0.13);
	EXPECT_DOUBLE_EQ(agent1.cell_adhesion_affinities()[3], 0.14);

	EXPECT_DOUBLE_EQ(agent0.chemotactic_sensitivities()[0], 2.0);
	EXPECT_DOUBLE_EQ(agent0.chemotactic_sensitivities()[1], 3.0);
	EXPECT_DOUBLE_EQ(agent1.chemotactic_sensitivities()[0], 4.0);
	EXPECT_DOUBLE_EQ(agent1.chemotactic_sensitivities()[1], 5.0);

	ASSERT_EQ(agent0.neighbors().size(), 1U);
	EXPECT_EQ(agent0.neighbors()[0], static_cast<index_t>(1));

	ASSERT_EQ(agent1.neighbors().size(), 2U);
	EXPECT_EQ(agent1.neighbors()[0], static_cast<index_t>(0));
	EXPECT_EQ(agent1.neighbors()[1], static_cast<index_t>(2));
}

} // namespace physicore::mechanics::physicell::tests
