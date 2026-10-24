#include "container_solver.h"

#include <algorithm>

#include "physicell/mechanical_agent_interface.h"

namespace physicore::mechanics::physicell::kernels::openmp_solver {

namespace {
void remove_attached(index_t to_remove, std::vector<index_t>* PHYSICORE_RESTRICT attached)
{
	for (auto spring : attached[to_remove])
	{
		auto it = std::ranges::find(attached[spring], to_remove);

		assert(it != attached[spring].end());

		*it = attached[spring].back();
		attached[spring].pop_back();
	}
}

void rename_attached(index_t old_index, index_t new_index, std::vector<index_t>* PHYSICORE_RESTRICT attached)
{
	for (auto spring : attached[old_index])
	{
		auto it = std::ranges::find(attached[spring], old_index);

		assert(it != attached[spring].end());

		*it = new_index;
	}
}

void remove_single(index_t i, const cell_flag_t* PHYSICORE_RESTRICT flag,
				   std::vector<index_t>* PHYSICORE_RESTRICT springs, mechanical_container_ptr& container,
				   index_t& counter)
{
	if (flag[i] != cell_flag_t::REMOVE)
		return;

	counter++;

	remove_attached(i, springs);
	rename_attached(container->size() - 1, i, springs);

	container->remove_at(i);
}

} // namespace

void container_solver::update_cell_container_for_phenotype(environment& e)
{
	auto& data = retrieve_agent_data(*e.agents);

#pragma omp master
	{
		for (std::size_t i = e.agents->size(); i > 0; i--)
		{
			remove_single(i - 1, data.state_data.flags.data(), data.state_data.springs.data(), e.agents,
						  e.removals_count);
		}
	}
#pragma omp barrier
}

} // namespace physicore::mechanics::physicell::kernels::openmp_solver
