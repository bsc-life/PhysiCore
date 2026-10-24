#include <biofvm/microenvironment.h>
#include <biofvm/microenvironment_builder.h>
#include <gtest/gtest.h>

using namespace physicore;
using namespace physicore::reactions_diffusion::biofvm;

namespace {
// Builds a 3x3x3 mesh with a single, non-diffusing, non-decaying substrate and pins the
// x-row at (x, 1, 1) to known values so the gradient can be computed analytically.
std::unique_ptr<microenvironment> gradient_test_microenv()
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 0.0, 0.0, 5.0);
	builder.resize(3, { 0, 0, 0 }, { 30, 30, 30 }, { 10, 10, 10 });

	builder.add_dirichlet_node({ 0, 1, 1 }, { 0.0 }, { true });
	builder.add_dirichlet_node({ 1, 1, 1 }, { 10.0 }, { true });
	builder.add_dirichlet_node({ 2, 1, 1 }, { 30.0 }, { true });

	auto env = builder.build();
	env->solver->initialize(*env);
	env->run_single_timestep();

	return env;
}
} // namespace

TEST(SubstrateGradientTest, CentralDifferenceAtInteriorVoxel)
{
	auto env = gradient_test_microenv();

	const auto gradient = env->solver->get_substrate_gradient(*env, 0, 1, 1, 1);

	EXPECT_NEAR(gradient[0], 1.5, 1e-9); // (30 - 0) / (2 * 10)
	EXPECT_NEAR(gradient[1], 0.0, 1e-9);
	EXPECT_NEAR(gradient[2], 0.0, 1e-9);
}

TEST(SubstrateGradientTest, OneSidedDifferenceAtMinBoundary)
{
	auto env = gradient_test_microenv();

	const auto gradient = env->solver->get_substrate_gradient(*env, 0, 0, 1, 1);

	EXPECT_NEAR(gradient[0], 1.0, 1e-9); // (10 - 0) / 10
}

TEST(SubstrateGradientTest, OneSidedDifferenceAtMaxBoundary)
{
	auto env = gradient_test_microenv();

	const auto gradient = env->solver->get_substrate_gradient(*env, 0, 2, 1, 1);

	EXPECT_NEAR(gradient[0], 2.0, 1e-9); // (30 - 10) / 10
}
