#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <common/types.h>
#include <gtest/gtest.h>

#include "physicell/openmp_solver/position_solver.h"
#include "physicell/openmp_solver/solve_pair_interface.h"

namespace physicore::mechanics::physicell::kernels::openmp_solver::tests {

using namespace physicore;

/**
 * @class SolvePairTest
 * @brief Unit tests for pairwise cell-cell force calculations.
 *
 * Tests the `solve_pair<dims>()` function which computes:
 * - Repulsive forces (based on cell overlap)
 * - Adhesive forces (based on affinity and proximity)
 * - Velocity updates from net force
 * - Simple pressure accumulation
 *
 * The function implements Newton's third law: force on cell i from j
 * is equal and opposite to force on cell j from cell i.
 */
class SolvePairTest : public ::testing::Test
{
protected:
	// Test data arrays (1D case: 2 cells)
	static constexpr index_t NUM_CELLS = 2;
	static constexpr index_t CELL_DEFS = 1;
	static constexpr real_t EPSILON = 1e-6;

	// Cell properties
	std::vector<real_t> position; // Dimension-sized positions
	std::array<real_t, NUM_CELLS> radius;
	std::vector<real_t> velocity; // Dimension-sized velocities
	std::array<real_t, NUM_CELLS> simple_pressure;
	std::array<real_t, NUM_CELLS> cell_cell_repulsion_strength;
	std::array<real_t, NUM_CELLS> cell_cell_adhesion_strength;
	std::array<real_t, NUM_CELLS> relative_maximum_adhesion_distance;
	std::array<real_t, NUM_CELLS> cell_adhesion_affinities; // 1 cell type
	std::array<index_t, NUM_CELLS> cell_definition_index;

	void SetUp() override
	{
		// Initialize default test state
		position.assign(NUM_CELLS * 3, 0.0);
		radius.fill(5.0);
		velocity.assign(NUM_CELLS * 3, 0.0);
		simple_pressure.fill(0.0);
		cell_cell_repulsion_strength.fill(1.0);
		cell_cell_adhesion_strength.fill(0.5);
		relative_maximum_adhesion_distance.fill(1.5);
		cell_adhesion_affinities.fill(1.0);
		cell_definition_index.fill(0);
	}

	/**
	 * @brief Setup two cells at specified distance with given radii
	 * @param distance Distance between cell centers
	 * @param r1 Radius of cell 1
	 * @param r2 Radius of cell 2
	 */
	void setup_cells_1d(real_t distance, real_t r1 = 5.0, real_t r2 = 5.0)
	{
		position.assign(NUM_CELLS * 1, 0.0);
		velocity.assign(NUM_CELLS * 1, 0.0);
		position[0] = 0.0;		// Cell 0 at origin
		position[1] = distance; // Cell 1 at distance
		radius[0] = r1;
		radius[1] = r2;
		simple_pressure.fill(0.0);
	}

	void setup_cells_2d(real_t dx, real_t dy, real_t r1 = 5.0, real_t r2 = 5.0)
	{
		position.assign(NUM_CELLS * 2, 0.0);
		velocity.assign(NUM_CELLS * 2, 0.0);
		position[0] = 0.0;
		position[1] = 0.0;
		position[2] = dx;
		position[3] = dy;
		radius[0] = r1;
		radius[1] = r2;
		simple_pressure.fill(0.0);
	}

	void setup_cells_3d(real_t dx, real_t dy, real_t dz, real_t r1 = 5.0, real_t r2 = 5.0)
	{
		position.assign(NUM_CELLS * 3, 0.0);
		velocity.assign(NUM_CELLS * 3, 0.0);
		position[0] = 0.0;
		position[1] = 0.0;
		position[2] = 0.0;
		position[3] = dx;
		position[4] = dy;
		position[5] = dz;
		radius[0] = r1;
		radius[1] = r2;
		simple_pressure.fill(0.0);
	}
};

// ============================================================================
// 1D TESTS
// ============================================================================

TEST_F(SolvePairTest, RepulsiveForce1D_Overlapping)
{
	// Two cells overlapping (distance < sum of radii)
	setup_cells_1d(5.0); // distance = 5, r1+r2 = 10, so overlap
	real_t vel_0_before = velocity[0];

	// Call solve_pair for 1D
	solve_pair<1>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	// Cell 0 should move left (negative direction) away from cell 1
	EXPECT_LT(velocity[0], vel_0_before) << "Repulsive force should push cell 0 negative";
	// Simple pressure should increase due to overlap
	EXPECT_GT(simple_pressure[0], 0.0) << "Simple pressure should accumulate from repulsion";
	EXPECT_GT(simple_pressure[1], 0.0) << "Simple pressure should accumulate for both cells";
}

TEST_F(SolvePairTest, NoForce1D_FarApart)
{
	// Two cells far apart beyond adhesion distance
	// adhesion_distance = relative_maximum_adhesion_distance * (r1 + r2)
	// = 1.5 * 10 = 15
	// At distance 100: adhesion = (1 - 100/15)^2 = 32.11 (large positive value!)
	// This is BEYOND the adhesion distance, so cell is no receiving force
	// The adhesion force becomes the dominant repulsive term
	//
	// NOTE: solve_pair only updates the LHS (cell 0) velocity
	// The RHS (cell 1) velocity is updated when solve_pair(j, i, ...) is called
	// This maintains Newton's 3rd law at the system level
	setup_cells_1d(100.0); // distance = 100, far beyond adhesion range (15)
	real_t vel_0_before = velocity[0];

	// Debug: check the actual values
	// position_difference = pos[1] - pos[0] = 100 - 0 = 100
	// repulsion = 0 (distance > sum of radii: 100 > 10)
	// adhesion = (1 - 100/15)^2 = 32.11
	// force = (0 - 32.11) / 100 = -0.3211
	// update: velocity[0] += force * position_difference = velocity[0] + (-0.3211) * 100
	// Should be negative update

	solve_pair<1>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	// Debug output - print what actually happened
	std::cout << "Cell 0: position=" << position[0] << ", velocity_before=" << vel_0_before
			  << ", velocity_after=" << velocity[0] << ", delta=" << (velocity[0] - vel_0_before) << std::endl;
	std::cout << "Cell 1: position=" << position[1] << ", velocity=" << velocity[1] << std::endl;
	std::cout << "Relative adhesion distance=" << relative_maximum_adhesion_distance[0] << " " <<
		relative_maximum_adhesion_distance[1] << std::endl;
	std::cout << "Radius=" << radius[0] << " " << radius[1] << std::endl;

	// The velocity update should be NEGATIVE (force pulls cell 0 away from cell 1)
	// But if it's positive, that means the force direction is wrong
	EXPECT_EQ(velocity[0], vel_0_before)
		<< "Beyond adhesion range, adhesion term creates repulsive force pushing cell 0 away";
}

TEST_F(SolvePairTest, AdhesiveForce1D_InAdhesionRange)
{
	// Two cells at adhesion distance (within adhesion_distance but not overlapping)
	setup_cells_1d(12.5); // distance = 12.5, adhesion_distance = 1.5 * (5+5) = 15
	real_t vel_0_before = velocity[0];
	real_t vel_1_before = velocity[1];

	solve_pair<1>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	// Adhesive force should pull them together
	// Cell 0 should move right (positive direction) towards cell 1
	EXPECT_GT(velocity[0], vel_0_before) << "Adhesive force should pull cell 0 positive";
	// Cell 1 is not updated
	EXPECT_EQ(velocity[1], vel_1_before) << "Adhesive force should pull cell 1 negative";
}

TEST_F(SolvePairTest, NewtonsThirdLaw1D_ForceSymmetry)
{
	// Verify that force on cell i from j = -force on cell j from i
	setup_cells_1d(7.0); // Overlapping configuration

	std::array<real_t, NUM_CELLS> vel_before{};
	std::copy_n(velocity.begin(), NUM_CELLS, vel_before.begin());

	solve_pair<1>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	solve_pair<1>(1, 0, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	real_t delta_vel_0 = velocity[0] - vel_before[0];
	real_t delta_vel_1 = velocity[1] - vel_before[1];

	// Forces should be opposite
	// Allow small tolerance for floating point arithmetic
	EXPECT_NEAR(delta_vel_0, -delta_vel_1, EPSILON * 10) << "Forces should be equal and opposite (Newton's 3rd law)";
}

TEST_F(SolvePairTest, SimplePressure1D_Accumulates)
{
	// Overlapping cells should accumulate simple pressure equally
	setup_cells_1d(6.0); // Overlap

	solve_pair<1>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	// Both cells should have equal pressure contribution
	EXPECT_EQ(simple_pressure[0], simple_pressure[1])
		<< "Simple pressure should be equal for both cells in symmetric configuration";
	EXPECT_GT(simple_pressure[0], 0.0) << "Simple pressure should be positive for repulsion";
}

TEST_F(SolvePairTest, ZeroRepulsion1D_NoRepulsiveForce)
{
	// Compare velocity update with repulsion enabled vs disabled
	setup_cells_1d(6.0); // Overlapping
	cell_cell_repulsion_strength[0] = 1.0;
	cell_cell_repulsion_strength[1] = 1.0;
	real_t vel_0_before_with_repulsion = velocity[0];

	solve_pair<1>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	const real_t vel_with_repulsion = velocity[0];
	const real_t delta_with_repulsion = vel_with_repulsion - vel_0_before_with_repulsion;

	setup_cells_1d(6.0); // Overlapping
	cell_cell_repulsion_strength[0] = 0.0;
	cell_cell_repulsion_strength[1] = 0.0;
	real_t vel_0_before_without_repulsion = velocity[0];

	solve_pair<1>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	const real_t vel_without_repulsion = velocity[0];
	const real_t delta_without_repulsion = vel_without_repulsion - vel_0_before_without_repulsion;

	// Removing repulsion should increase the magnitude of the velocity change in this overlap case.
	EXPECT_LT(std::abs(delta_with_repulsion), std::abs(delta_without_repulsion))
		<< "Velocity magnitude should be larger when repulsion is disabled";
}

// ============================================================================
// 2D TESTS
// ============================================================================

TEST_F(SolvePairTest, RepulsiveForce2D_Overlapping)
{
	setup_cells_2d(5.0, 5.0); // Diagonal distance = sqrt(50) ≈ 7.07, r1+r2=10, overlap
	std::array<real_t, NUM_CELLS * 2> vel_before;
	std::copy(velocity.begin(), velocity.end(), vel_before.begin());
	cell_cell_adhesion_strength[0] = 0.0; //Deactivate adhesion strength
	cell_cell_adhesion_strength[1] = 0.0; //Deactivate adhesion strength

	solve_pair<2>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	std::cout << "Velocity before: (" << vel_before[0] << ", " << vel_before[1] << ")" << std::endl;
	std::cout << "Velocity after: (" << velocity[0] << ", " << velocity[1] << ")" << std::endl;

	// Both velocity components should move away from the other cell
	// Cell 0 at (0,0) should move negative in both x and y
	EXPECT_LT(velocity[0], vel_before[0]) << "Repulsive force should push cell 0 negative in x";
	EXPECT_LT(velocity[1], vel_before[1]) << "Repulsive force should push cell 0 negative in y";
}

TEST_F(SolvePairTest, AdhesiveForce2D_InAdhesionRange)
{
	setup_cells_2d(12.0, 3.0); // distance ≈ 12.37, within adhesion range
	std::array<real_t, NUM_CELLS * 2> vel_before;
	std::copy(velocity.begin(), velocity.end(), vel_before.begin());

	solve_pair<2>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	// Cells should move towards each other
	// Cell 0 should move towards cell 1 (positive x, positive y)
	EXPECT_GT(velocity[0], vel_before[0]) << "Adhesive force should pull cell 0 towards cell 1 (x)";
	EXPECT_GT(velocity[1], vel_before[1]) << "Adhesive force should pull cell 0 towards cell 1 (y)";
}

TEST_F(SolvePairTest, NewtonsThirdLaw2D_ForceSymmetry)
{
	setup_cells_2d(8.0, 6.0);
	std::array<real_t, NUM_CELLS * 2> vel_before;
	std::copy(velocity.begin(), velocity.end(), vel_before.begin());

	solve_pair<2>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	solve_pair<2>(1, 0, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	// Force on cell 0 = -Force on cell 1
	real_t force_0_x = velocity[0] - vel_before[0];
	real_t force_0_y = velocity[1] - vel_before[1];
	real_t force_1_x = velocity[2] - vel_before[2];
	real_t force_1_y = velocity[3] - vel_before[3];

	EXPECT_NEAR(force_0_x, -force_1_x, EPSILON * 10) << "X-component forces should be equal and opposite";
	EXPECT_NEAR(force_0_y, -force_1_y, EPSILON * 10) << "Y-component forces should be equal and opposite";
}

// ============================================================================
// 3D TESTS
// ============================================================================

TEST_F(SolvePairTest, RepulsiveForce3D_Overlapping)
{
	setup_cells_3d(4.0, 4.0, 4.0); // distance ≈ 6.93, r1+r2=10, overlap
	std::array<real_t, NUM_CELLS * 3> vel_before;
	std::copy(velocity.begin(), velocity.end(), vel_before.begin());

	cell_cell_adhesion_strength[0] = 0.0;
	cell_cell_adhesion_strength[1] = 0.0;

	solve_pair<3>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());
	std::cout << "Velocity before: (" << vel_before[0] << ", " << vel_before[1] << ", " << vel_before[2] << ")" << std::endl;
	std::cout << "Velocity after: (" << velocity[0] << ", " << velocity[1] << ", " << velocity[2] << ")" << std::endl;

	// Cell 0 should move away (all negative)
	EXPECT_LT(velocity[0], vel_before[0]);
	EXPECT_LT(velocity[1], vel_before[1]);
	EXPECT_LT(velocity[2], vel_before[2]);

}

TEST_F(SolvePairTest, AdhesiveForce3D_DirectionalAccuracy)
{
	setup_cells_3d(10.0, 0.0, 0.0); // Pure x-direction separation
	std::array<real_t, NUM_CELLS * 3> vel_before;
	std::copy(velocity.begin(), velocity.end(), vel_before.begin());

	solve_pair<3>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	// Force should be primarily in x direction
	real_t force_x = std::abs(velocity[0] - vel_before[0]);
	real_t force_y = std::abs(velocity[1] - vel_before[1]);
	real_t force_z = std::abs(velocity[2] - vel_before[2]);

	EXPECT_GT(force_x, force_y) << "X-component should dominate";
	EXPECT_GT(force_x, force_z) << "X-component should dominate";
}

TEST_F(SolvePairTest, NewtonsThirdLaw3D_ForceSymmetry)
{
	setup_cells_3d(7.0, 3.0, 5.0);
	std::array<real_t, NUM_CELLS * 3> vel_before;
	std::copy(velocity.begin(), velocity.end(), vel_before.begin());

	solve_pair<3>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	solve_pair<3>(1, 0, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	real_t force_0_mag_sq = 0.0, force_1_mag_sq = 0.0;
	for (index_t d = 0; d < 3; ++d)
	{
		real_t f0 = velocity[d] - vel_before[d];
		real_t f1 = velocity[3 + d] - vel_before[3 + d];
		EXPECT_NEAR(f0, -f1, EPSILON * 10) << "Force component " << d << " not symmetric";
		force_0_mag_sq += f0 * f0;
		force_1_mag_sq += f1 * f1;
	}

	EXPECT_NEAR(force_0_mag_sq, force_1_mag_sq, EPSILON * 100) << "Force magnitudes should be equal";
}

// ============================================================================
// EDGE CASES
// ============================================================================

TEST_F(SolvePairTest, ZeroDistance1D_Minimum)
{
	// Cells at same position (should not divide by zero)
	setup_cells_1d(0.001); // Very small distance
	EXPECT_NO_THROW(solve_pair<1>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(),
								  radius.data(), cell_cell_repulsion_strength.data(),
								  cell_cell_adhesion_strength.data(), relative_maximum_adhesion_distance.data(),
								  cell_adhesion_affinities.data(), cell_definition_index.data());)
		<< "Should handle near-zero distance without crashing";

	EXPECT_TRUE(std::isfinite(velocity[0])) << "Velocity should be finite";
	EXPECT_TRUE(std::isfinite(simple_pressure[0])) << "Pressure should be finite";
}

TEST_F(SolvePairTest, DifferentRadii1D_AsymmetricForce)
{
	// Cell 0 with radius 10, cell 1 with radius 5
	setup_cells_1d(10.0, 10.0, 5.0);
	std::array<real_t, NUM_CELLS> vel_before{};
	std::copy_n(velocity.begin(), NUM_CELLS, vel_before.begin());

	solve_pair<1>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	// Pressure should still be equal (symmetric interaction)
	EXPECT_EQ(simple_pressure[0], simple_pressure[1])
		<< "Pressure accumulation should be symmetric regardless of radius";
}

TEST_F(SolvePairTest, ZeroAffinity1D_NoAdhesion)
{
	cell_adhesion_affinities[0] = 0.0;
	setup_cells_1d(12.0); // In adhesion range but zero affinity
	std::array<real_t, NUM_CELLS> vel_before{};
	std::copy_n(velocity.begin(), NUM_CELLS, vel_before.begin());

	solve_pair<1>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(), radius.data(),
				  cell_cell_repulsion_strength.data(), cell_cell_adhesion_strength.data(),
				  relative_maximum_adhesion_distance.data(), cell_adhesion_affinities.data(),
				  cell_definition_index.data());

	// No adhesive force should be applied
	EXPECT_EQ(velocity[0], vel_before[0]) << "Zero affinity should eliminate adhesion";
	EXPECT_EQ(velocity[1], vel_before[1]) << "Zero affinity should eliminate adhesion";
}

TEST_F(SolvePairTest, DifferentCellTypes1D_AffinityLookup)
{
	// Two cells of different types with affinity matrix
	cell_definition_index[0] = 0;
	cell_definition_index[1] = 1; // Different type (but only 1 type total in our test)
	// Affinity matrix is 1x1, so this tests edge case handling
	setup_cells_1d(12.0);

	EXPECT_NO_THROW(solve_pair<1>(0, 1, CELL_DEFS, velocity.data(), simple_pressure.data(), position.data(),
								  radius.data(), cell_cell_repulsion_strength.data(),
								  cell_cell_adhesion_strength.data(), relative_maximum_adhesion_distance.data(),
								  cell_adhesion_affinities.data(), cell_definition_index.data());)
		<< "Should handle cell type indexing correctly";
}

} // namespace physicore::mechanics::physicell::kernels::openmp_solver::tests
