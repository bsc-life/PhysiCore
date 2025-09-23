#include "solver_registry.h"

#include <cassert>

#include "cpu_solver/register_solver.h"

using namespace physicore::biofvm;

void solver_registry::register_factory(std::string solver_name, solver_factory_func_t&& f)
{
	if (factory_registry.contains(solver_name))
	{
		assert(false);
		return;
	}

	factory_registry.emplace(std::move(solver_name), std::move(f));
}

std::unique_ptr<solver> solver_registry::get(const std::string& solver_name)
{
	if (!factory_registry.contains(solver_name))
	{
		assert(false);
		return nullptr;
	}

	return factory_registry[solver_name]();
}

solver_registry& solver_registry::instance()
{
	static solver_registry r;
	return r;
}

struct attachment_point
{
	attachment_point() { kernels::cpu::attach_to_registry(); }
};

static attachment_point ap;
