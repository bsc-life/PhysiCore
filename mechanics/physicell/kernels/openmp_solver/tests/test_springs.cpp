#include <gtest/gtest.h>

#include "position_solver.h"

using namespace physicore;
using namespace physicore::mechanics::physicell;

class UpdateSpringAttachmentsComplexTest : public ::testing::TestWithParam<index_t>
{};

class agent_retriever : public physicore::generic_agent_solver<mechanical_agent>
{
public:
	using physicore::generic_agent_solver<mechanical_agent>::retrieve_agent_data;
};

TEST(UpdateSpringAttachmentsTest, Simple2D)
{
	const index_t dims = 2;
	environment env({ dims, { -500, -500, -500 }, { 500, 500, 500 }, { 20, 20, 20 } }, 1, 1, 0.1);

	auto create_agent = [&](real_t x, real_t y) {
		auto* agent = env.agents->create();
		agent->radius() = 9;
		agent->is_movable() = 1;
		agent->attachment_elastic_constant() = 0.01;
		agent->detachment_rate() = 0;
		agent->cell_adhesion_affinities()[0] = 2;
		agent->position()[0] = x;
		agent->position()[1] = y;
		return agent;
	};
	auto* a1 = create_agent(0, 0);
	auto* a2 = create_agent(0, 100);
	auto* a3 = create_agent(100, 0);

	auto& data = agent_retriever().retrieve_agent_data(*env.agents);
	data.state_data.springs[0] = { 1, 2 };
	data.state_data.springs[1] = { 0, 2 };
	data.state_data.springs[2] = { 0, 1 };

	kernels::openmp_solver::position_solver solver;
	solver.update_spring_attachments(env);
	solver.update_positions(env);

	EXPECT_FLOAT_EQ(a1->position()[0], 0.3);
	EXPECT_FLOAT_EQ(a1->position()[1], 0.3);

	EXPECT_FLOAT_EQ(a2->position()[0], 0.3);
	EXPECT_FLOAT_EQ(a2->position()[1], 99.4);

	EXPECT_FLOAT_EQ(a3->position()[0], 99.4);
	EXPECT_FLOAT_EQ(a3->position()[1], 0.3);

	solver.update_spring_attachments(env);
	solver.update_positions(env);

	EXPECT_FLOAT_EQ(a1->position()[0], 0.4973);
	EXPECT_FLOAT_EQ(a1->position()[1], 0.4973);

	EXPECT_FLOAT_EQ(a2->position()[0], 0.4973);
	EXPECT_FLOAT_EQ(a2->position()[1], 99.0054);

	EXPECT_FLOAT_EQ(a3->position()[0], 99.0054);
	EXPECT_FLOAT_EQ(a3->position()[1], 0.4973);
}

TEST(UpdateSpringAttachmentsTest, Complex2D)
{
	const index_t dims = 2;
	environment env({ dims, { -500, -500, -500 }, { 500, 500, 500 }, { 20, 20, 20 } }, 2, 1, 0.1);

	auto create_agent = [&](real_t x, real_t y, index_t type) {
		auto* agent = env.agents->create();
		agent->radius() = 9;
		agent->is_movable() = 1;
		agent->attachment_elastic_constant() = type == 0 ? 0.01 : 0.02;
		agent->detachment_rate() = 0;
		agent->cell_adhesion_affinities()[0] = (real_t)type * 2 + 2;
		agent->cell_adhesion_affinities()[1] = (real_t)type * 2 + 3;
		agent->position()[0] = x;
		agent->position()[1] = y;
		agent->agent_type_index() = type;
		return agent;
	};
	auto* a1 = create_agent(0, 0, 0);
	auto* a2 = create_agent(0, 100, 0);
	auto* a3 = create_agent(100, 0, 1);

	auto& data = agent_retriever().retrieve_agent_data(*env.agents);
	data.state_data.springs[0] = { 1, 2 };
	data.state_data.springs[1] = { 0, 2 };
	data.state_data.springs[2] = { 0, 1 };

	kernels::openmp_solver::position_solver solver;
	solver.update_spring_attachments(env);
	solver.update_positions(env);

	EXPECT_FLOAT_EQ(a1->position()[0], 0.734847);
	EXPECT_FLOAT_EQ(a1->position()[1], 0.3);

	EXPECT_FLOAT_EQ(a2->position()[0], 0.734847);
	EXPECT_FLOAT_EQ(a2->position()[1], 98.96515);

	EXPECT_FLOAT_EQ(a3->position()[0], 98.5303);
	EXPECT_FLOAT_EQ(a3->position()[1], 0.734847);
}

TEST(UpdateSpringAttachmentsTest, Simple3D)
{
	const index_t dims = 3;
	environment env({ dims, { -500, -500, -500 }, { 500, 500, 500 }, { 20, 20, 20 } }, 1, 1, 0.1);

	auto create_agent = [&](real_t x, real_t y, real_t z) {
		auto* agent = env.agents->create();
		agent->radius() = 9;
		agent->is_movable() = 1;
		agent->attachment_elastic_constant() = 0.01;
		agent->detachment_rate() = 0;
		agent->cell_adhesion_affinities()[0] = 2;
		agent->position()[0] = x;
		agent->position()[1] = y;
		agent->position()[2] = z;
		return agent;
	};
	auto* a1 = create_agent(0, 0, 0);
	auto* a2 = create_agent(0, 100, 0);
	auto* a3 = create_agent(100, 0, 0);
	auto* a4 = create_agent(0, 0, 100);

	auto& data = agent_retriever().retrieve_agent_data(*env.agents);
	data.state_data.springs[0] = { 1, 2, 3 };
	data.state_data.springs[1] = { 0, 2, 3 };
	data.state_data.springs[2] = { 0, 1, 3 };
	data.state_data.springs[3] = { 0, 1, 2 };

	kernels::openmp_solver::position_solver solver;
	solver.update_spring_attachments(env);
	solver.update_positions(env);

	EXPECT_FLOAT_EQ(a1->position()[0], 0.3);
	EXPECT_FLOAT_EQ(a1->position()[1], 0.3);
	EXPECT_FLOAT_EQ(a1->position()[2], 0.3);

	EXPECT_FLOAT_EQ(a2->position()[0], 0.3);
	EXPECT_FLOAT_EQ(a2->position()[1], 99.1);
	EXPECT_FLOAT_EQ(a2->position()[2], 0.3);

	EXPECT_FLOAT_EQ(a3->position()[0], 99.1);
	EXPECT_FLOAT_EQ(a3->position()[1], 0.3);
	EXPECT_FLOAT_EQ(a3->position()[2], 0.3);

	EXPECT_FLOAT_EQ(a4->position()[0], 0.3);
	EXPECT_FLOAT_EQ(a4->position()[1], 0.3);
	EXPECT_FLOAT_EQ(a4->position()[2], 99.1);

	solver.update_spring_attachments(env);
	solver.update_positions(env);

	EXPECT_FLOAT_EQ(a1->position()[0], 0.4964);
	EXPECT_FLOAT_EQ(a1->position()[1], 0.4964);
	EXPECT_FLOAT_EQ(a1->position()[2], 0.4964);

	EXPECT_FLOAT_EQ(a2->position()[0], 0.4964);
	EXPECT_FLOAT_EQ(a2->position()[1], 98.5108);
	EXPECT_FLOAT_EQ(a2->position()[2], 0.4964);

	EXPECT_FLOAT_EQ(a3->position()[0], 98.5108);
	EXPECT_FLOAT_EQ(a3->position()[1], 0.4964);
	EXPECT_FLOAT_EQ(a3->position()[2], 0.4964);

	EXPECT_FLOAT_EQ(a4->position()[0], 0.4964);
	EXPECT_FLOAT_EQ(a4->position()[1], 0.4964);
	EXPECT_FLOAT_EQ(a4->position()[2], 98.5108);
}

TEST_P(UpdateSpringAttachmentsComplexTest, NoMove)
{
	const index_t dims = GetParam();
	environment env({ dims, { -500, -500, -500 }, { 500, 500, 500 }, { 20, 20, 20 } }, 1, 1, 0.1);

	auto create_agent = [&](real_t x) {
		auto* agent = env.agents->create();
		agent->radius() = 9;
		agent->is_movable() = 1;
		agent->attachment_elastic_constant() = 0.01;
		agent->detachment_rate() = 0;
		agent->cell_adhesion_affinities()[0] = 2;
		for (index_t d = 0; d < dims; ++d)
			agent->position()[d] = x;
		return agent;
	};

	auto* a1 = create_agent(0);
	auto* a2 = create_agent(100);

	auto& data = agent_retriever().retrieve_agent_data(*env.agents);
	data.state_data.springs[0] = { 1 };
	data.state_data.springs[1] = { 0 };

	env.automated_spring_adhesion = false;

	kernels::openmp_solver::position_solver solver;
	solver.update_spring_attachments(env);
	solver.update_positions(env);

	for (index_t d = 0; d < dims; ++d)
	{
		EXPECT_FLOAT_EQ(a1->position()[d], 0);
		EXPECT_FLOAT_EQ(a2->position()[d], 100);
	}

	env.automated_spring_adhesion = true;
	a1->is_movable() = 0;
	a2->is_movable() = 0;
	solver.update_spring_attachments(env);
	solver.update_positions(env);

	for (index_t d = 0; d < dims; ++d)
	{
		EXPECT_FLOAT_EQ(a1->position()[d], 0);
		EXPECT_FLOAT_EQ(a2->position()[d], 100);
	}
}

TEST_P(UpdateSpringAttachmentsComplexTest, AttachAndDetach)
{
	const index_t dims = GetParam();
	environment env({ dims, { -500, -500, -500 }, { 500, 500, 500 }, { 20, 20, 20 } }, 1, 1, 0.1);

	auto create_agent = [&](real_t x) {
		auto* agent = env.agents->create();
		agent->radius() = 9;
		agent->is_movable() = 1;
		agent->attachment_elastic_constant() = 0.01;
		agent->attachment_rate() = 0;
		agent->detachment_rate() = 0;
		agent->cell_adhesion_affinities()[0] = 2;
		agent->maximum_number_of_attachments() = 1;
		for (index_t d = 0; d < dims; ++d)
			agent->position()[d] = x;
		return agent;
	};

	auto* a1 = create_agent(0);
	auto* a2 = create_agent(100);

	auto& data = agent_retriever().retrieve_agent_data(*env.agents);
	data.state_data.neighbors[0] = { 1 };
	data.state_data.neighbors[1] = { 0 };

	a1->attachment_rate() = 1000;

	kernels::openmp_solver::position_solver solver;
	solver.update_spring_attachments(env);
	solver.update_positions(env);

	// Just one cell has attachment rate set, but both cells should attach
	ASSERT_EQ(data.state_data.springs[0].size(), 1);
	ASSERT_EQ(data.state_data.springs[1].size(), 1);
	EXPECT_EQ(data.state_data.springs[0][0], 1);
	EXPECT_EQ(data.state_data.springs[1][0], 0);

	for (index_t d = 0; d < dims; ++d)
	{
		EXPECT_FLOAT_EQ(a1->position()[d], 0.3);
		EXPECT_FLOAT_EQ(a2->position()[d], 99.7);
	}

	a1->attachment_rate() = 0;
	a2->detachment_rate() = 1000;

	solver.update_spring_attachments(env);
	solver.update_positions(env);

	EXPECT_EQ(data.state_data.springs[0].size(), 0);
	EXPECT_EQ(data.state_data.springs[1].size(), 0);

	for (index_t d = 0; d < dims; ++d)
	{
		EXPECT_FLOAT_EQ(a1->position()[d], 0.2);
		EXPECT_FLOAT_EQ(a2->position()[d], 99.8);
	}
}

TEST_P(UpdateSpringAttachmentsComplexTest, MaxAttachmentsLimit)
{
	const index_t dims = GetParam();
	environment env({ dims, { -500, -500, -500 }, { 500, 500, 500 }, { 20, 20, 20 } }, 1, 1, 0.1);

	auto create_agent = [&](real_t x) {
		auto* agent = env.agents->create();
		agent->radius() = 9;
		agent->is_movable() = 1;
		agent->attachment_elastic_constant() = 0.01;
		agent->attachment_rate() = 1000;
		agent->detachment_rate() = 0;
		agent->cell_adhesion_affinities()[0] = 2;
		for (index_t d = 0; d < dims; ++d)
			agent->position()[d] = x;
		return agent;
	};

	auto* a1 = create_agent(0);
	auto* a2 = create_agent(100);
	auto* a3 = create_agent(200);

	a1->maximum_number_of_attachments() = 2;
	a2->maximum_number_of_attachments() = 2;
	a3->maximum_number_of_attachments() = 1;

	auto& data = agent_retriever().retrieve_agent_data(*env.agents);
	data.state_data.neighbors[0] = { 1, 2 };
	data.state_data.neighbors[1] = { 0, 2 };
	data.state_data.neighbors[2] = { 0, 1 };

	kernels::openmp_solver::position_solver solver;
	solver.update_spring_attachments(env);
	solver.update_positions(env);

	ASSERT_EQ(data.state_data.springs[2].size(), 1);
	ASSERT_GE(data.state_data.springs[0].size(), 1);
	ASSERT_GE(data.state_data.springs[1].size(), 1);
	ASSERT_EQ(data.state_data.springs[0].size() + data.state_data.springs[1].size(), 3);
}

INSTANTIATE_TEST_SUITE_P(Dimensions, UpdateSpringAttachmentsComplexTest, ::testing::Values(index_t(2), index_t(3)));
