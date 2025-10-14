#include <gtest/gtest.h>

#include "solver.h"
#include "solver_registry.h"

using namespace physicore;
using namespace physicore::biofvm;

class mock_solver : public solver
{
public:
	void solve(microenvironment&, index_t) override {}
	real_t get_substrate_density(index_t, index_t, index_t, index_t) const override { return 0; }
};

TEST(SolverRegistryTest, CheckPresentSolvers)
{
	auto registry = solver_registry::instance();

	EXPECT_NE(registry.get("openmp_solver"), nullptr);

#ifdef PHYSICORE_HAS_THRUST
	EXPECT_NE(registry.get("thrust_solver"), nullptr);
#endif
}

TEST(SolverRegistryTest, GetAndSet)
{
	solver_registry registry;

	// Add a new one
	{
		EXPECT_EQ(registry.register_factory("solver_x", []() { return std::make_unique<mock_solver>(); }), true);

		EXPECT_NE(registry.get("solver_x"), nullptr);
	}

	// Try to add the same one again
	{
		EXPECT_EQ(registry.register_factory("solver_x", []() { return std::make_unique<mock_solver>(); }), false);
	}

	// Test non-existing
	{
#if NDEBUG or _DEBUG
		EXPECT_EQ(registry.get("solver_y"), nullptr);
#endif
	}
}

TEST(SolverRegistryTest, RegistryAdder)
{
	auto& registry = solver_registry::instance();

	registry_adder<mock_solver> openmp_solver_adder("solver_x");

	EXPECT_NE(registry.get("solver_x"), nullptr);
}
