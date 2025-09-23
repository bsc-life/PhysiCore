#include <array>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

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

TEST(MicroenvironmentBuilder, BulkFunctionsAndInternalizedSubstrates)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	builder.resize(3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 });

	bool bulk_called = false;
	builder.set_bulk_functions(
		[&](index_t, index_t, index_t, index_t) {
			bulk_called = true;
			return 0;
		},
		nullptr, nullptr);
	builder.do_compute_internalized_substrates();

	auto env = builder.build();
	ASSERT_TRUE(env->compute_internalized_substrates);
	ASSERT_TRUE(env->supply_rate_func != nullptr);

	// Call bulk function to check assignment
	env->supply_rate_func(0, 0, 0, 0);
	ASSERT_TRUE(bulk_called);
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
}
