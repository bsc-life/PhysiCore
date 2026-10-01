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
				   index_t& counter, std::mutex& m)
{
	while (true)
	{
		if (flag[i] != cell_flag_t::REMOVE)
			return;

		{
			const std::unique_lock<std::mutex> l(m);

			// this can happen, it is safe because we are not deallocating memory for removed elements
			if (i >= container->size())
				return;

			if (flag[i] == cell_flag_t::REMOVE)
				counter++;

			remove_attached(i, springs);
			rename_attached(container->size() - 1, i, springs);

			container->remove_at(i);

			if (i == container->size())
				return;
		}
	}
}

void divide_single(index_t i, cell_flag_t* PHYSICORE_RESTRICT flag, const index_t* PHYSICORE_RESTRICT agent_type_index,
				   environment& e, index_t& counter, std::mutex& m)
{
	if (flag[i] == cell_flag_t::DIVIDE)
	{
		flag[i] = cell_flag_t::NONE;

		{
			const std::unique_lock<std::mutex> l(m);

			e.create_with_type(agent_type_index[i]);
			counter++;
		}
	}
}

} // namespace

void container_solver::update_cell_container_for_phenotype(environment& e)
{
	auto& data = retrieve_agent_data(*e.agents);

	const auto n = e.agents->size();
#pragma omp barrier

#pragma omp for
	for (std::size_t i = 0; i < n; i++)
	{
		divide_single(i, data.state_data.flags.data(), data.state_data.agent_type_index.data(), e, e.divisions_count,
					  container_mutex);
	}

#pragma omp for
	for (std::size_t i = 0; i < n; i++)
	{
		remove_single(i, data.state_data.flags.data(), data.state_data.springs.data(), e.agents, e.removals_count,
					  container_mutex);
	}
}

} // namespace physicore::mechanics::physicell::kernels::openmp_solver
