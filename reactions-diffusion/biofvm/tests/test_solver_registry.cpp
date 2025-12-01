#include <gtest/gtest.h>

#include "solver.h"
#include "solver_registry.h"

using namespace physicore;
using namespace physicore::biofvm;

class mock_solver : public solver
{
public:
	void initialize(microenvironment& /*m*/) override {}
	void solve(microenvironment& /*m*/, index_t /*iterations*/) override {}
	real_t get_substrate_density(index_t /*s*/, index_t /*x*/, index_t /*y*/, index_t /*z*/) const override
	{
		return 0;
	}
	real_t& get_substrate_density(index_t /*s*/, index_t /*x*/, index_t /*y*/, index_t /*z*/) override
	{
		static real_t dummy = 0;
		return dummy;
	}
	void reinitialize_dirichlet(microenvironment& /*m*/) override {}
};

TEST(SolverRegistryTest, CheckPresentSolvers)
{
	auto registry = solver_registry::instance();

	EXPECT_NE(registry.get("openmp_solver"), nullptr);

#ifdef PHYSICORE_HAS_TBB_THRUST
	EXPECT_NE(registry.get("tbb_thrust_solver"), nullptr);
#endif
#ifdef PHYSICORE_HAS_CUDA_THRUST
	EXPECT_NE(registry.get("cuda_thrust_solver"), nullptr);
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
#ifdef NDEBUG
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
