#include <gtest/gtest.h>

#include "micromechanics/solver.h"
#include "micromechanics/solver_registry.h"

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

class mock_solver : public solver
{
public:
	void initialize(environment& /*e*/) override {}
	void update_cell_neighbors(environment& /*e*/) override {}
	void update_cell_forces(environment& /*e*/) override {}
	void calculate_cell_data(environment& /*e*/) override {}
	void update_motility(environment& /*e*/) override {}
	void update_basement_membrane_interactions(environment& /*e*/) override {}
	void update_spring_attachments(environment& /*e*/) override {}
	void update_positions(environment& /*e*/) override {}
};

TEST(SolverRegistryTest, CheckPresentSolvers)
{
	auto& registry = solver_registry::instance();

	// The singleton registry should be usable (register + construct a solver).
	registry_adder<mock_solver> const mock_solver_adder("test_mock_solver_present");
	EXPECT_TRUE(registry.is_available("test_mock_solver_present"));
	EXPECT_NE(registry.get("test_mock_solver_present"), nullptr);
}

TEST(SolverRegistryTest, GetAndSet)
{
	solver_registry registry;

	// Add a new solver
	{
		EXPECT_TRUE(registry.register_factory("mock_solver_a", []() { return std::make_unique<mock_solver>(); }));
		EXPECT_TRUE(registry.is_available("mock_solver_a"));
		EXPECT_NE(registry.get("mock_solver_a"), nullptr);
	}

	// Try to add the same one again - should fail
	{
		EXPECT_FALSE(registry.register_factory("mock_solver_a", []() { return std::make_unique<mock_solver>(); }));
	}

	// Test non-existing solver
	{
#ifdef NDEBUG
		EXPECT_EQ(registry.get("nonexistent_solver"), nullptr);
		EXPECT_FALSE(registry.is_available("nonexistent_solver"));
#endif
	}
}

TEST(SolverRegistryTest, RegistryAdder)
{
	auto& registry = solver_registry::instance();

	// Use registry_adder to register a solver
	registry_adder<mock_solver> const mock_solver_adder("test_mock_solver");

	EXPECT_TRUE(registry.is_available("test_mock_solver"));
	EXPECT_NE(registry.get("test_mock_solver"), nullptr);
}

TEST(SolverRegistryTest, AvailableSolvers)
{
	solver_registry registry;

	EXPECT_TRUE(registry.available_solvers().empty());

	registry.register_factory("solver_1", []() { return std::make_unique<mock_solver>(); });
	registry.register_factory("solver_2", []() { return std::make_unique<mock_solver>(); });

	auto solvers = registry.available_solvers();
	EXPECT_EQ(solvers.size(), 2U);

	// Check both solvers are in the list (order may vary due to unordered_map)
	bool found_1 = false;
	bool found_2 = false;
	for (const auto& name : solvers)
	{
		if (name == "solver_1")
			found_1 = true;
		if (name == "solver_2")
			found_2 = true;
	}
	EXPECT_TRUE(found_1);
	EXPECT_TRUE(found_2);
}

TEST(SolverRegistryTest, CreatedSolversAreUnique)
{
	solver_registry registry;
	registry.register_factory("unique_solver", []() { return std::make_unique<mock_solver>(); });

	auto solver1 = registry.get("unique_solver");
	auto solver2 = registry.get("unique_solver");

	EXPECT_NE(solver1, nullptr);
	EXPECT_NE(solver2, nullptr);
	EXPECT_NE(solver1.get(), solver2.get()); // Different instances
}
