#include <biofvm/microenvironment_builder.h>
#include <gtest/gtest.h>

using namespace physicore;
using namespace physicore::biofvm;

TEST(DirichletModification, SetInteriorVoxel)
{
	// Build a simple microenvironment with one interior Dirichlet voxel
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	builder.add_density("Glucose", "mM", 0.5, 0.02, 5.0);
	builder.resize(3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 });

	// Add a Dirichlet node at (5, 5, 5)
	std::vector<real_t> values = { 100.0, 50.0 };
	std::vector<bool> conditions = { true, true };
	builder.add_dirichlet_node({ 5, 5, 5 }, values, conditions);

	auto env = builder.build();

	// Verify initial values
	EXPECT_EQ(env->dirichlet_interior_voxels_count, 1);
	EXPECT_EQ(env->dirichlet_interior_values[0], 100.0); // O2 value
	EXPECT_EQ(env->dirichlet_interior_values[1], 50.0);	 // Glucose value
	EXPECT_TRUE(env->dirichlet_interior_conditions[0]);	 // O2 enabled
	EXPECT_TRUE(env->dirichlet_interior_conditions[1]);	 // Glucose enabled

	// Modify the O2 value at voxel 0
	EXPECT_NO_THROW(env->update_dirichlet_interior_voxel({ 5, 5, 5 }, 0, 200.0, true));
	EXPECT_EQ(env->dirichlet_interior_values[0], 200.0);
	EXPECT_TRUE(env->dirichlet_interior_conditions[0]);

	// Modify the Glucose value and disable it
	EXPECT_NO_THROW(env->update_dirichlet_interior_voxel({ 5, 5, 5 }, 1, 75.0, false));
	EXPECT_EQ(env->dirichlet_interior_values[1], 75.0);
	EXPECT_FALSE(env->dirichlet_interior_conditions[1]);
}

TEST(DirichletModification, SetInteriorVoxelOutOfBounds)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	builder.resize(3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 });

	// Add one Dirichlet node
	builder.add_dirichlet_node({ 5, 5, 5 }, { 100.0 }, { true });

	auto env = builder.build();

	// Try to set substrate index out of bounds
	EXPECT_THROW(env->update_dirichlet_interior_voxel({ 5, 5, 5 }, 1, 200.0, true), std::runtime_error);
}

TEST(DirichletModification, SetBoundaryMin)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	builder.add_density("Glucose", "mM", 0.5, 0.02, 5.0);
	builder.resize(3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 });

	// Add boundary Dirichlet conditions via builder
	builder.add_boundary_dirichlet_conditions(0,					   // O2
											  { 100.0, 110.0, 120.0 }, // min values for x, y, z
											  { 0.0, 0.0, 0.0 },	   // max values (not used here)
											  { true, true, true },	   // min conditions enabled
											  { false, false, false }  // max conditions disabled
	);

	auto env = builder.build();

	// Verify initial values
	EXPECT_EQ(env->dirichlet_min_boundary_values[0][0], 100.0); // x-min O2
	EXPECT_EQ(env->dirichlet_min_boundary_values[1][0], 110.0); // y-min O2
	EXPECT_EQ(env->dirichlet_min_boundary_values[2][0], 120.0); // z-min O2
	EXPECT_TRUE(env->dirichlet_min_boundary_conditions[0][0]);
	EXPECT_TRUE(env->dirichlet_min_boundary_conditions[1][0]);
	EXPECT_TRUE(env->dirichlet_min_boundary_conditions[2][0]);

	// Modify x-min boundary for O2
	EXPECT_NO_THROW(env->update_dirichlet_boundary_min('x', 0, 150.0, true));
	EXPECT_EQ(env->dirichlet_min_boundary_values[0][0], 150.0);
	EXPECT_TRUE(env->dirichlet_min_boundary_conditions[0][0]);

	// Modify y-min boundary for Glucose (which wasn't set initially)
	EXPECT_NO_THROW(env->update_dirichlet_boundary_min('y', 1, 25.0, true));
	EXPECT_EQ(env->dirichlet_min_boundary_values[1][1], 25.0);
	EXPECT_TRUE(env->dirichlet_min_boundary_conditions[1][1]);

	// Disable z-min boundary for O2
	EXPECT_NO_THROW(env->update_dirichlet_boundary_min('z', 0, 120.0, false));
	EXPECT_EQ(env->dirichlet_min_boundary_values[2][0], 120.0);
	EXPECT_FALSE(env->dirichlet_min_boundary_conditions[2][0]);
}

TEST(DirichletModification, SetBoundaryMax)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	builder.resize(3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 });

	// Add boundary Dirichlet conditions via builder
	builder.add_boundary_dirichlet_conditions(0,					   // O2
											  { 0.0, 0.0, 0.0 },	   // min values (not used here)
											  { 200.0, 210.0, 220.0 }, // max values for x, y, z
											  { false, false, false }, // min conditions disabled
											  { true, true, true }	   // max conditions enabled
	);

	auto env = builder.build();

	// Verify initial values
	EXPECT_EQ(env->dirichlet_max_boundary_values[0][0], 200.0);
	EXPECT_EQ(env->dirichlet_max_boundary_values[1][0], 210.0);
	EXPECT_EQ(env->dirichlet_max_boundary_values[2][0], 220.0);
	EXPECT_TRUE(env->dirichlet_max_boundary_conditions[0][0]);
	EXPECT_TRUE(env->dirichlet_max_boundary_conditions[1][0]);
	EXPECT_TRUE(env->dirichlet_max_boundary_conditions[2][0]);

	// Modify x-max boundary
	EXPECT_NO_THROW(env->update_dirichlet_boundary_max('x', 0, 250.0, true));
	EXPECT_EQ(env->dirichlet_max_boundary_values[0][0], 250.0);
	EXPECT_TRUE(env->dirichlet_max_boundary_conditions[0][0]);

	// Disable y-max boundary
	EXPECT_NO_THROW(env->update_dirichlet_boundary_max('y', 0, 210.0, false));
	EXPECT_EQ(env->dirichlet_max_boundary_values[1][0], 210.0);
	EXPECT_FALSE(env->dirichlet_max_boundary_conditions[1][0]);
}

TEST(DirichletModification, SetBoundaryOutOfBounds)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	builder.resize(2, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 }); // 2D mesh

	auto env = builder.build();

	// Try to set dimension out of bounds (z-axis in 2D mesh)
	EXPECT_THROW(env->update_dirichlet_boundary_min('z', 0, 100.0, true), std::runtime_error);
	EXPECT_THROW(env->update_dirichlet_boundary_max('z', 0, 100.0, true), std::runtime_error);

	// Try to set substrate index out of bounds
	EXPECT_THROW(env->update_dirichlet_boundary_min('x', 1, 100.0, true), std::runtime_error);
	EXPECT_THROW(env->update_dirichlet_boundary_max('x', 1, 100.0, true), std::runtime_error);
}

TEST(DirichletModification, EmptyInteriorVoxel)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	builder.add_density("Glucose", "mM", 0.5, 0.02, 5.0);
	builder.resize(3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 });

	// Build without any boundary conditions
	auto env = builder.build();

	EXPECT_EQ(env->dirichlet_interior_voxels_count, 0);

	// Set an interior voxel that doesn't exist
	env->update_dirichlet_interior_voxel({ 5, 5, 5 }, 0, 100.0, true);

	// Verify the value we set
	EXPECT_EQ(env->dirichlet_interior_voxels[0], 5);
	EXPECT_EQ(env->dirichlet_interior_voxels[1], 5);
	EXPECT_EQ(env->dirichlet_interior_voxels[2], 5);
	EXPECT_EQ(env->dirichlet_interior_values[0], 100.0);
	EXPECT_TRUE(env->dirichlet_interior_conditions[0]);
	// Verify other substrates are initialized to defaults
	EXPECT_EQ(env->dirichlet_interior_values[1], 0.0);
	EXPECT_FALSE(env->dirichlet_interior_conditions[1]);
}

TEST(DirichletModification, LazyBoundaryAllocation)
{
	microenvironment_builder builder;
	builder.add_density("O2", "mmHg", 1.0, 0.01, 20.0);
	builder.add_density("Glucose", "mM", 0.5, 0.02, 5.0);
	builder.resize(3, { 0, 0, 0 }, { 10, 10, 10 }, { 1, 1, 1 });

	// Build without any boundary conditions
	auto env = builder.build();

	// Verify that boundary arrays are initially null
	EXPECT_EQ(env->dirichlet_min_boundary_values[0], nullptr);
	EXPECT_EQ(env->dirichlet_max_boundary_values[0], nullptr);

	// Set a boundary condition - this should allocate the arrays
	EXPECT_NO_THROW(env->update_dirichlet_boundary_min('x', 0, 100.0, true));

	// Verify arrays are now allocated
	EXPECT_NE(env->dirichlet_min_boundary_values[0], nullptr);
	EXPECT_NE(env->dirichlet_min_boundary_conditions[0], nullptr);

	// Verify the value we set
	EXPECT_EQ(env->dirichlet_min_boundary_values[0][0], 100.0);
	EXPECT_TRUE(env->dirichlet_min_boundary_conditions[0][0]);

	// Verify other substrates are initialized to defaults
	EXPECT_EQ(env->dirichlet_min_boundary_values[0][1], 0.0);
	EXPECT_FALSE(env->dirichlet_min_boundary_conditions[0][1]);
}
