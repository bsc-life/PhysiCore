#include <array>
#include <string>

#include <gtest/gtest.h>

#include "position_solver.h"
#include "reactions_diffusion/reactions_diffusion_interface.h"

using namespace physicore;
using namespace physicore::mechanics::physicell;

class mock_diffusion_interface : public reactions_diffusion::reactions_diffusion_interface
{
public:
	std::array<std::array<real_t, 3>, 2> gradients {};

	std::span<const std::string> get_substrate_names() const override { return {}; }
	std::span<const std::string> get_substrate_units() const override { return {}; }

	real_t get_substrate_density(index_t /*s*/, std::span<const real_t> /*position*/) const override { return 0; }

	std::array<real_t, 3> get_substrate_gradient(index_t substrate, std::span<const real_t> /*position*/) const override
	{
		return gradients[substrate];
	}

	void run_single_timestep() override {}
	void serialize_state(real_t /*current_time*/) override {}
};

namespace {
mechanical_agent_interface* setup_chemotaxis_environment(environment& env,
														 std::shared_ptr<mock_diffusion_interface>& diffusion)
{
	diffusion = std::make_shared<mock_diffusion_interface>();
	env.diffusion = diffusion;
	return env.agents->create();
}
} // namespace

class agent_retriever : public generic_agent_solver<mechanical_agent>
{
public:
	using generic_agent_solver<mechanical_agent>::retrieve_agent_data;
};

class UpdateMotilityTest : public ::testing::TestWithParam<index_t>
{};

TEST_P(UpdateMotilityTest, NonMotileAgentIsUnchanged)
{
	const index_t dims = GetParam();
	environment env({ dims, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 } }, 1, 1, 0);
	auto* agent = env.agents->create();

	for (index_t d = 0; d < dims; ++d)
	{
		agent->motility_vector()[d] = static_cast<real_t>(d + 1);
		agent->velocity()[d] = static_cast<real_t>(10 + d);
	}

	kernels::openmp_solver::position_solver solver;
	solver.update_motility(env);

	for (index_t d = 0; d < dims; ++d)
	{
		EXPECT_FLOAT_EQ(agent->motility_vector()[d], static_cast<real_t>(d + 1));
		EXPECT_FLOAT_EQ(agent->velocity()[d], static_cast<real_t>(10 + d));
	}
}

TEST(MigrationBiasFunctorTest, NoneReturnsNull)
{
	environment env({ 2, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 } }, 1, 2, 0);
	std::shared_ptr<mock_diffusion_interface> diffusion;
	setup_chemotaxis_environment(env, diffusion);

	kernels::openmp_solver::position_solver solver;

	EXPECT_EQ(solver.create_migration_bias_functor(env, migration_bias_type::NONE), nullptr);
}

TEST(UpdateMotilityTest, SimpleChemotaxis2D)
{
	environment env({ 2, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 } }, 1, 2, 0.1);
	std::shared_ptr<mock_diffusion_interface> diffusion;
	auto* agent = setup_chemotaxis_environment(env, diffusion);
	diffusion->gradients[1] = { 1, 2, 3 };
	agent->chemotaxis_index() = 1;
	agent->chemotaxis_direction() = 1;
	agent->migration_speed() = 4;
	agent->persistence_time() = 0;
	agent->migration_bias() = 1;
	agent->is_motile() = 1;

	kernels::openmp_solver::position_solver solver;
	agent->migration_bias_functor() = solver.create_migration_bias_functor(env, migration_bias_type::SIMPLE);

	for (index_t iters = 0; iters < 4; iters++)
	{
		solver.update_motility(env);
		EXPECT_FLOAT_EQ(agent->motility_vector()[0], 1.7888544);
		EXPECT_FLOAT_EQ(agent->motility_vector()[1], 3.5777087);
	}
}

TEST(UpdateMotilityTest, SimpleChemotaxis3D)
{
	environment env({ 3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 } }, 1, 2, 0.1);
	std::shared_ptr<mock_diffusion_interface> diffusion;
	auto* agent = setup_chemotaxis_environment(env, diffusion);
	diffusion->gradients[1] = { 1, 2, 3 };
	agent->chemotaxis_index() = 1;
	agent->chemotaxis_direction() = -1;
	agent->migration_speed() = 4;
	agent->persistence_time() = 0;
	agent->migration_bias() = 1;
	agent->is_motile() = 1;

	kernels::openmp_solver::position_solver solver;
	agent->migration_bias_functor() = solver.create_migration_bias_functor(env, migration_bias_type::SIMPLE);

	for (index_t iters = 0; iters < 4; iters++)
	{
		solver.update_motility(env);
		EXPECT_FLOAT_EQ(agent->motility_vector()[0], -1.0690449);
		EXPECT_FLOAT_EQ(agent->motility_vector()[1], -2.13809);
		EXPECT_FLOAT_EQ(agent->motility_vector()[2], -3.207135);
	}
}

TEST(UpdateMotilityTest, AdvancedChemotaxis2D)
{
	environment env({ 2, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 } }, 1, 2, 0.1);
	std::shared_ptr<mock_diffusion_interface> diffusion;
	auto* agent = setup_chemotaxis_environment(env, diffusion);
	diffusion->gradients[0] = { 1, 2, 3 };
	diffusion->gradients[1] = { 4, 5, 6 };
	agent->chemotactic_sensitivities()[0] = 7;
	agent->chemotactic_sensitivities()[1] = -8;
	agent->chemotaxis_direction() = 1;
	agent->migration_speed() = 4;
	agent->persistence_time() = 0;
	agent->migration_bias() = 1;
	agent->is_motile() = 1;

	kernels::openmp_solver::position_solver solver;
	agent->migration_bias_functor() = solver.create_migration_bias_functor(env, migration_bias_type::ADVANCED);

	solver.update_motility(env);
	EXPECT_FLOAT_EQ(agent->motility_vector()[0], -2.7724349);
	EXPECT_FLOAT_EQ(agent->motility_vector()[1], -2.8833323);

	diffusion->gradients[0] = { 7, 8, 9 };
	diffusion->gradients[1] = { 10, 11, 12 };

	solver.update_motility(env);
	EXPECT_FLOAT_EQ(agent->motility_vector()[0], -2.78318);
	EXPECT_FLOAT_EQ(agent->motility_vector()[1], -2.87296);
}

TEST(UpdateMotilityTest, AdvancedChemotaxis3D)
{
	environment env({ 3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 } }, 1, 2, 0.1);
	std::shared_ptr<mock_diffusion_interface> diffusion;
	auto* agent = setup_chemotaxis_environment(env, diffusion);
	diffusion->gradients[0] = { 1, 2, 3 };
	diffusion->gradients[1] = { 4, 5, 6 };
	agent->chemotactic_sensitivities()[0] = 7;
	agent->chemotactic_sensitivities()[1] = -8;
	agent->chemotaxis_direction() = 1;
	agent->migration_speed() = 4;
	agent->persistence_time() = 0;
	agent->migration_bias() = 1;
	agent->is_motile() = 1;

	kernels::openmp_solver::position_solver solver;
	agent->migration_bias_functor() = solver.create_migration_bias_functor(env, migration_bias_type::ADVANCED);

	solver.update_motility(env);
	EXPECT_FLOAT_EQ(agent->motility_vector()[0], -2.2194839);
	EXPECT_FLOAT_EQ(agent->motility_vector()[1], -2.3082631);
	EXPECT_FLOAT_EQ(agent->motility_vector()[2], -2.3970425);

	diffusion->gradients[0] = { 7, 8, 9 };
	diffusion->gradients[1] = { 10, 11, 12 };

	solver.update_motility(env);
	EXPECT_FLOAT_EQ(agent->motility_vector()[0], -2.2365043);
	EXPECT_FLOAT_EQ(agent->motility_vector()[1], -2.30865);
	EXPECT_FLOAT_EQ(agent->motility_vector()[2], -2.380795);
}

TEST(UpdateMotilityTest, AdvancedChemotaxisNormalized2D)
{
	environment env({ 2, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 } }, 1, 2, 0.1);
	std::shared_ptr<mock_diffusion_interface> diffusion;
	auto* agent = setup_chemotaxis_environment(env, diffusion);
	diffusion->gradients[0] = { 1, 2, 3 };
	diffusion->gradients[1] = { 4, 5, 6 };
	agent->chemotactic_sensitivities()[0] = 7;
	agent->chemotactic_sensitivities()[1] = -8;
	agent->chemotaxis_direction() = 1;
	agent->migration_speed() = 4;
	agent->persistence_time() = 0;
	agent->migration_bias() = 1;
	agent->is_motile() = 1;

	kernels::openmp_solver::position_solver solver;
	agent->migration_bias_functor() =
		solver.create_migration_bias_functor(env, migration_bias_type::ADVANCED_NORMALIZED);

	solver.update_motility(env);
	EXPECT_FLOAT_EQ(agent->motility_vector()[0], -3.999887);
	EXPECT_FLOAT_EQ(agent->motility_vector()[1], 0.030078145);

	diffusion->gradients[0] = { 7, 8, 9 };
	diffusion->gradients[1] = { 10, 11, 12 };

	solver.update_motility(env);
	EXPECT_FLOAT_EQ(agent->motility_vector()[0], -3.0567069);
	EXPECT_FLOAT_EQ(agent->motility_vector()[1], -2.5800278);
}

TEST(UpdateMotilityTest, AdvancedChemotaxisNormalized3D)
{
	environment env({ 3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 } }, 1, 2, 0.1);
	std::shared_ptr<mock_diffusion_interface> diffusion;
	auto* agent = setup_chemotaxis_environment(env, diffusion);
	diffusion->gradients[0] = { 1, 2, 3 };
	diffusion->gradients[1] = { 4, 5, 6 };
	agent->chemotactic_sensitivities()[0] = 7;
	agent->chemotactic_sensitivities()[1] = -8;
	agent->chemotaxis_direction() = 1;
	agent->migration_speed() = 4;
	agent->persistence_time() = 0;
	agent->migration_bias() = 1;
	agent->is_motile() = 1;

	kernels::openmp_solver::position_solver solver;
	agent->migration_bias_functor() =
		solver.create_migration_bias_functor(env, migration_bias_type::ADVANCED_NORMALIZED);

	solver.update_motility(env);
	EXPECT_FLOAT_EQ(agent->motility_vector()[0], -3.6244786);
	EXPECT_FLOAT_EQ(agent->motility_vector()[1], -1.6669483);
	EXPECT_FLOAT_EQ(agent->motility_vector()[2], 0.290582);

	diffusion->gradients[0] = { 7, 8, 9 };
	diffusion->gradients[1] = { 10, 11, 12 };

	solver.update_motility(env);
	EXPECT_FLOAT_EQ(agent->motility_vector()[0], -2.62217);
	EXPECT_FLOAT_EQ(agent->motility_vector()[1], -2.2937832);
	EXPECT_FLOAT_EQ(agent->motility_vector()[2], -1.965397);
}

INSTANTIATE_TEST_SUITE_P(Dimensions, UpdateMotilityTest, ::testing::Values(index_t(1), index_t(2), index_t(3)));
