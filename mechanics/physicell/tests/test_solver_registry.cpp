#include <memory>
#include <string>

#include <common/types.h>
#include <gtest/gtest.h>

#include "physicell/solver_registry.h"

using namespace physicore;
using namespace physicore::mechanics::physicell;

namespace {

class DummySolver final : public solver
{
public:
	void initialize(environment& /*e*/) override {}
	void solve(environment& /*e*/, index_t /*iterations*/) override {}
};

} // namespace

TEST(SolverRegistryTest, InstanceReturnsSameRegistry)
{
	auto& r1 = solver_registry::instance();
	auto& r2 = solver_registry::instance();
	EXPECT_EQ(&r1, &r2);
}

TEST(SolverRegistryTest, RegisterFactoryStoresAndCreatesSolver)
{
	auto& registry = solver_registry::instance();
	const std::string name = "test_solver_registry_dummy";

	const bool registered = registry.register_factory(name, []() { return std::make_unique<DummySolver>(); });
	EXPECT_TRUE(registered);

	auto solver_instance = registry.get(name);
	ASSERT_NE(solver_instance, nullptr);
	EXPECT_NE(dynamic_cast<DummySolver*>(solver_instance.get()), nullptr);
}

TEST(SolverRegistryTest, RegisterFactoryRejectsDuplicate)
{
	auto& registry = solver_registry::instance();
	const std::string name = "test_solver_registry_duplicate";

	const bool first = registry.register_factory(name, []() { return std::make_unique<DummySolver>(); });
	const bool second = registry.register_factory(name, []() { return std::make_unique<DummySolver>(); });

	EXPECT_TRUE(first);
	EXPECT_FALSE(second);

	auto solver_instance = registry.get(name);
	ASSERT_NE(solver_instance, nullptr);
}

#ifdef NDEBUG
TEST(SolverRegistryTest, GetReturnsNullptrForMissingKeyInRelease)
{
	auto& registry = solver_registry::instance();
	EXPECT_EQ(registry.get("test_solver_registry_missing"), nullptr);
}
#endif
