#include <gtest/gtest.h>
#include <physicell/openmp_solver/position_solver.h>

#include "common/cartesian_mesh.h"
#include "physicell/environment.h"

using namespace physicore;
using namespace physicore::mechanics::physicell;

class UpdateBasementMembraneTest : public ::testing::TestWithParam<index_t>
{};

TEST_P(UpdateBasementMembraneTest, SimpleEdge)
{
	const index_t dims = GetParam();
	environment env(0.1, dims, 1, 1);
	env.set_mesh(physicore::cartesian_mesh { dims, { -500, -500, -500 }, { 500, 500, 500 }, { 20, 20, 20 } });

	auto* a1 = env.agents->create();
	a1->cell_BM_repulsion_strength() = 100;
	a1->radius() = 9;
	a1->is_movable() = 1;

	for (index_t d = 0; d < dims; ++d)
		a1->position()[d] = -500;

	kernels::openmp_solver::position_solver solver;
	solver.update_basement_membrane_interactions(env, env.get_mesh());
	solver.update_positions(env);

	for (index_t d = 0; d < dims; ++d)
		EXPECT_FLOAT_EQ(a1->position()[d], -485);

	for (index_t i = 0; i < 10; ++i)
	{
		solver.update_basement_membrane_interactions(env, env.get_mesh());
		solver.update_positions(env);

		for (index_t d = 0; d < dims; ++d)
			EXPECT_FLOAT_EQ(a1->position()[d], -490);
	}
}

TEST_P(UpdateBasementMembraneTest, MultipleEdge)
{
	const index_t dims = GetParam();
	environment env(0.1, dims, 1, 1);
	env.set_mesh(physicore::cartesian_mesh { dims, { -500, -500, -500 }, { 500, 500, 500 }, { 20, 20, 20 } });

	auto* a1 = env.agents->create();
	a1->cell_BM_repulsion_strength() = 100;
	a1->radius() = 9;
	a1->is_movable() = 1;

	auto* a2 = env.agents->create();
	a2->cell_BM_repulsion_strength() = 100;
	a2->radius() = 9;
	a2->is_movable() = 1;

	for (index_t d = 0; d < dims; ++d)
	{
		a1->position()[d] = -500;
		a2->position()[d] = 500;
	}

	kernels::openmp_solver::position_solver solver;
	solver.update_basement_membrane_interactions(env, env.get_mesh());
	solver.update_positions(env);

	for (index_t d = 0; d < dims; ++d)
	{
		EXPECT_FLOAT_EQ(a1->position()[d], -485);
		EXPECT_FLOAT_EQ(a2->position()[d], 485);
	}

	for (index_t i = 0; i < 10; ++i)
	{
		solver.update_basement_membrane_interactions(env, env.get_mesh());
		solver.update_positions(env);

		for (index_t d = 0; d < dims; ++d)
		{
			EXPECT_FLOAT_EQ(a1->position()[d], -490);
			EXPECT_FLOAT_EQ(a2->position()[d], 490);
		}
	}
}

TEST_P(UpdateBasementMembraneTest, SimpleCenter)
{
	const index_t dims = GetParam();
	environment env(0.1, dims, 1, 1);
	env.set_mesh(physicore::cartesian_mesh { dims, { -500, -500, -500 }, { 500, 500, 500 }, { 20, 20, 20 } });

	auto* a1 = env.agents->create();
	a1->cell_BM_repulsion_strength() = 100;
	a1->radius() = 9;
	a1->is_movable() = 1;

	for (index_t d = 0; d < dims; ++d)
		a1->position()[d] = 0;

	kernels::openmp_solver::position_solver solver;

	for (index_t i = 0; i < 10; ++i)
	{
		solver.update_basement_membrane_interactions(env, env.get_mesh());
		solver.update_positions(env);

		for (index_t d = 0; d < dims; ++d)
			EXPECT_FLOAT_EQ(a1->position()[d], 0);
	}
}

TEST_P(UpdateBasementMembraneTest, SimpleOneOff)
{
	const index_t dims = GetParam();
	environment env(0.1, dims, 1, 1);
	env.set_mesh(physicore::cartesian_mesh { dims, { -500, -500, -500 }, { 500, 500, 500 }, { 20, 20, 20 } });

	auto* a1 = env.agents->create();
	a1->cell_BM_repulsion_strength() = 100;
	a1->radius() = 9;
	a1->is_movable() = 1;

	for (index_t d = 0; d < dims; ++d)
		a1->position()[d] = d == dims - 1 ? 0 : -500;

	kernels::openmp_solver::position_solver solver;
	solver.update_basement_membrane_interactions(env, env.get_mesh());
	solver.update_positions(env);

	for (index_t d = 0; d < dims; ++d)
		EXPECT_FLOAT_EQ(a1->position()[d], d == dims - 1 ? 0 : -485);

	for (index_t i = 0; i < 10; ++i)
	{
		solver.update_basement_membrane_interactions(env, env.get_mesh());
		solver.update_positions(env);

		for (index_t d = 0; d < dims; ++d)
			EXPECT_FLOAT_EQ(a1->position()[d], d == dims - 1 ? 0 : -490);
	}
}

TEST_P(UpdateBasementMembraneTest, NoMove)
{
	const index_t dims = GetParam();
	environment env(0.1, dims, 1, 1);
	env.set_mesh(physicore::cartesian_mesh { dims, { -500, -500, -500 }, { 500, 500, 500 }, { 20, 20, 20 } });

	auto* a1 = env.agents->create();
	a1->cell_BM_repulsion_strength() = 100;
	a1->radius() = 9;
	a1->is_movable() = 1;

	for (index_t d = 0; d < dims; ++d)
		a1->position()[d] = -500;

	kernels::openmp_solver::position_solver solver;
	env.virtual_wall_at_domain_edges = false;
	solver.update_basement_membrane_interactions(env, env.get_mesh());
	solver.update_positions(env);

	for (index_t d = 0; d < dims; ++d)
		EXPECT_FLOAT_EQ(a1->position()[d], -500);

	env.virtual_wall_at_domain_edges = true;
	a1->is_movable() = 0;
	solver.update_basement_membrane_interactions(env, env.get_mesh());
	solver.update_positions(env);

	for (index_t d = 0; d < dims; ++d)
		EXPECT_FLOAT_EQ(a1->position()[d], -500);
}

INSTANTIATE_TEST_SUITE_P(Dimensions, UpdateBasementMembraneTest, ::testing::Values(index_t(1), index_t(2), index_t(3)));
