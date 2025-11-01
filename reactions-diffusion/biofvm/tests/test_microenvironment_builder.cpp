#include <array>
#include <memory>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "bulk_functor.h"
#include "microenvironment.h"
#include "microenvironment_builder.h"

using namespace physicore;
using namespace physicore::biofvm;

TEST(MicroenvironmentBuilder, SettersAndBuild)
{
	microenvironment_builder builder;
	builder.set_name("env");
	builder.set_time_units("min");
	builder.set_space_units("um");
	builder.set_time_step(0.1);

	builder.resize(3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 });
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	builder.add_density("Glucose", "mM", 0.5, 0.02, 5.0);

	auto env = builder.build();
	ASSERT_EQ(env->name, "env");
	ASSERT_EQ(env->time_units, "min");
	ASSERT_EQ(env->space_units, "um");
	ASSERT_EQ(env->substrates_names.size(), 2);
	EXPECT_EQ(env->substrates_names[0], "O2");
	EXPECT_EQ(env->substrates_names[1], "Glucose");
	EXPECT_NE(env->initial_conditions, nullptr);
	EXPECT_NE(env->diffusion_coefficients, nullptr);
	EXPECT_NE(env->decay_rates, nullptr);
}

TEST(MicroenvironmentBuilder, GetDensityIndex)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	builder.add_density("Glucose", "mM", 0.5, 0.02, 5.0);

	EXPECT_EQ(builder.get_density_index("O2"), 0);
	EXPECT_EQ(builder.get_density_index("Glucose"), 1);
	EXPECT_THROW(builder.get_density_index("NotFound"), std::runtime_error);
}

TEST(MicroenvironmentBuilder, AddDirichletNode)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	builder.resize(3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 });

	std::vector<real_t> values = { 1.0 };
	std::vector<bool> conditions = { true };
	// Should not throw
	EXPECT_NO_THROW(builder.add_dirichlet_node({ 0, 0, 0 }, values, conditions));

	// Build and verify dirichlet node is present
	auto env = builder.build();
	ASSERT_EQ(env->dirichlet_interior_voxels_count, 1);
	// Check that the first voxel index is correct (should be 0 for {0,0,0})
	EXPECT_EQ(env->dirichlet_interior_voxels[0], 0);
	// Check that the value and condition match what we set
	EXPECT_EQ(env->dirichlet_interior_values[0], values[0]);
	EXPECT_EQ(env->dirichlet_interior_conditions[0], conditions[0]);
}

TEST(MicroenvironmentBuilder, AddDirichletNodeImplicitCondition)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	builder.resize(3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 });

	std::vector<real_t> values = { 1.0 };
	std::vector<bool> conditions = { true };
	// Should not throw
	EXPECT_NO_THROW(builder.add_dirichlet_node({ 0, 0, 0 }, values));

	// Build and verify dirichlet node is present
	auto env = builder.build();
	ASSERT_EQ(env->dirichlet_interior_voxels_count, 1);
	// Check that the first voxel index is correct (should be 0 for {0,0,0})
	EXPECT_EQ(env->dirichlet_interior_voxels[0], 0);
	// Check that the value and condition match what we set
	EXPECT_EQ(env->dirichlet_interior_values[0], values[0]);
	EXPECT_EQ(env->dirichlet_interior_conditions[0], conditions[0]);
}

TEST(MicroenvironmentBuilder, AddDirichletNodeThrows)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);

	// No mesh set
	std::vector<real_t> values = { 1.0 };
	std::vector<bool> conditions = { true };
	EXPECT_THROW(builder.add_dirichlet_node({ 0, 0, 0 }, values, conditions), std::runtime_error);

	// Set mesh
	builder.resize(3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 });

	// Wrong values size
	std::vector<real_t> bad_values = { 1.0, 2.0 };
	EXPECT_THROW(builder.add_dirichlet_node({ 0, 0, 0 }, bad_values, conditions), std::runtime_error);

	// Wrong conditions size
	std::vector<bool> bad_conditions = { true, false };
	EXPECT_THROW(builder.add_dirichlet_node({ 0, 0, 0 }, values, bad_conditions), std::runtime_error);
}

TEST(MicroenvironmentBuilder, AddBoundaryDirichletConditions)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	builder.add_density("Glucose", "mM", 0.5, 0.02, 5.0);

	builder.resize(3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 });

	std::array<real_t, 3> mins_values_o2 = { 0.0, 1.0, 2.0 };
	std::array<real_t, 3> maxs_values_o2 = { 1.0, 2.0, 3.0 };
	std::array<bool, 3> mins_conditions_o2 = { true, false, false };
	std::array<bool, 3> maxs_conditions_02 = { false, true, false };

	builder.add_boundary_dirichlet_conditions(0, mins_values_o2, maxs_values_o2, mins_conditions_o2,
											  maxs_conditions_02);

	std::array<real_t, 3> mins_values_gl = { 10.0, 10.0, 10.0 };
	std::array<real_t, 3> maxs_values_gl = { 11.0, 11.0, 11.0 };
	std::array<bool, 3> mins_conditions_gl = { false, false, false };
	std::array<bool, 3> maxs_conditions_gl = { false, true, true };

	builder.add_boundary_dirichlet_conditions(1, mins_values_gl, maxs_values_gl, mins_conditions_gl,
											  maxs_conditions_gl);

	auto env = builder.build();

	ASSERT_THAT(std::span<real_t>(env->dirichlet_min_boundary_values[0].get(), 2), testing::ElementsAre(0.0, 10.0));
	ASSERT_EQ(env->dirichlet_min_boundary_values[1], nullptr);
	ASSERT_EQ(env->dirichlet_min_boundary_values[2], nullptr);

	ASSERT_EQ(env->dirichlet_max_boundary_values[0], nullptr);
	ASSERT_THAT(std::span<real_t>(env->dirichlet_max_boundary_values[1].get(), 2), testing::ElementsAre(2.0, 11.0));
	ASSERT_THAT(std::span<real_t>(env->dirichlet_max_boundary_values[2].get(), 2), testing::ElementsAre(3.0, 11.0));

	ASSERT_THAT(std::span<bool>(env->dirichlet_min_boundary_conditions[0].get(), 2), testing::ElementsAre(true, false));
	ASSERT_EQ(env->dirichlet_min_boundary_conditions[1], nullptr);
	ASSERT_EQ(env->dirichlet_min_boundary_conditions[2], nullptr);

	ASSERT_EQ(env->dirichlet_max_boundary_conditions[0], nullptr);
	ASSERT_THAT(std::span<bool>(env->dirichlet_max_boundary_conditions[1].get(), 2), testing::ElementsAre(true, true));
	ASSERT_THAT(std::span<bool>(env->dirichlet_max_boundary_conditions[2].get(), 2), testing::ElementsAre(false, true));
}

TEST(MicroenvironmentBuilder, AddBoundaryDirichletConditionsThrows)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);

	std::array<real_t, 3> mins_values = { 0.0, 0.0, 0.0 };
	std::array<real_t, 3> maxs_values = { 1.0, 1.0, 1.0 };
	std::array<bool, 3> mins_conditions = { true, false, false };
	std::array<bool, 3> maxs_conditions = { false, true, false };

	// Out of bounds index
	EXPECT_THROW(
		builder.add_boundary_dirichlet_conditions(1, mins_values, maxs_values, mins_conditions, maxs_conditions),
		std::runtime_error);

	// Valid index
	EXPECT_NO_THROW(
		builder.add_boundary_dirichlet_conditions(0, mins_values, maxs_values, mins_conditions, maxs_conditions));
}

struct test_functor : bulk_functor
{
	virtual real_t supply_rates(index_t, index_t, index_t, index_t) { return 42; }
	virtual real_t uptake_rates(index_t, index_t, index_t, index_t) { return 1; }
	virtual real_t supply_target_densities(index_t, index_t, index_t, index_t) { return 2; }
};

TEST(MicroenvironmentBuilder, BulkFunctionsAndInternalizedSubstrates)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	builder.resize(3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 });

	builder.set_bulk_functions(std::make_unique<test_functor>());
	builder.do_compute_internalized_substrates();

	auto env = builder.build();
	ASSERT_TRUE(env->compute_internalized_substrates);
	ASSERT_TRUE(env->bulk_fnc != nullptr);

	// Call bulk function to check assignment
	auto ret = env->bulk_fnc->supply_rates(0, 0, 0, 0);
	ASSERT_EQ(ret, 42);
}

TEST(MicroenvironmentBuilder, BuildThrows)
{
	microenvironment_builder builder;
	// No mesh
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	EXPECT_THROW(builder.build(), std::runtime_error);

	// No densities
	microenvironment_builder builder2;
	builder2.resize(3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 });
	EXPECT_THROW(builder2.build(), std::runtime_error);

#ifdef NDEBUG
	microenvironment_builder builder3;
	builder3.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	builder3.resize(3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 });
	builder3.select_solver("nonexistent");
	EXPECT_THROW(builder3.build(), std::runtime_error);
#endif
}
