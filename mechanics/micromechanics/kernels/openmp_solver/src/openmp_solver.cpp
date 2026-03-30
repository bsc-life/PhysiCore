#include "openmp_solver.h"

#include <tuple>

#include <common/base_agent_data.h>
#include <micromechanics/agent_container.h>
#include <micromechanics/agent_data.h>
#include <micromechanics/cell_aggregation.h>
#include <micromechanics/environment.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void openmp_solver::initialize(environment& e)
{
	if (initialized_)
		return;

	n_solver_.initialize(e);
	f_solver_.initialize(e);
	m_solver_.initialize(e);
	bm_solver_.initialize(e);
	s_solver_.initialize(e);
	p_solver_.initialize(e);

	initialized_ = true;
}

void openmp_solver::update_cell_neighbors(environment& e) { n_solver_.update_neighbors(e); }

void openmp_solver::update_cell_forces(environment& e) { f_solver_.calculate_forces(e); }

void openmp_solver::calculate_cell_data(environment& e)
{
	auto& agents = *e.agents;
	auto& base_data = *std::get<std::unique_ptr<base_agent_data>>(agents.agent_datas);
	auto& mech_data = *std::get<std::unique_ptr<agent_data>>(agents.agent_datas);

	aggregate_cell_data_from_agents(base_data, mech_data, e.cells);
}

void openmp_solver::update_motility(environment& e) { m_solver_.update_motility(e); }

void openmp_solver::update_basement_membrane_interactions(environment& e) { bm_solver_.update_interactions(e); }

void openmp_solver::update_spring_attachments(environment& e) { s_solver_.update_spring_attachments(e); }

void openmp_solver::update_positions(environment& e) { p_solver_.update_positions(e); }

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
