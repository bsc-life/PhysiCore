#include "openmp_solver.h"

#include <cmath>
#include <numbers>
#include <tuple>
#include <utility>

#include <common/base_agent_data.h>
#include <micromechanics/agent_data.h>
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

	const auto count = static_cast<index_t>(agents.size());
	constexpr index_t dims = 3;
	constexpr real_t four_thirds_pi = 4.0 / 3.0 * std::numbers::pi;

	// Clear previous cell data
	e.cells.clear();

	// Temporary accumulators for position and velocity calculation
	std::unordered_map<index_t, std::array<real_t, 3>> position_sums;
	std::unordered_map<index_t, index_t> position_counts;
	std::unordered_map<index_t, std::array<real_t, 3>> velocity_sums;
	std::unordered_map<index_t, index_t> velocity_counts;

	// First pass: accumulate all agent data per cell
	for (index_t i = 0; i < count; ++i)
	{
		index_t const cell_id = mech_data.cell_ids[i];

		// Skip standalone agents (cell_id == -1)
		if (std::cmp_equal(cell_id ,-1)))
			continue;

		std::uint8_t const agent_type = mech_data.agent_types[i];

		// Compartment counts
		e.cells.compartment_counts[{ cell_id, agent_type }]++;

		// Volume: sum of 4/3 * pi * r^3
		real_t const r = mech_data.radii[i];
		e.cells.volumes[cell_id] += four_thirds_pi * r * r * r;

		// Position accumulation
		real_t const px = base_data.positions[i * dims];
		real_t const py = base_data.positions[i * dims + 1];
		real_t const pz = base_data.positions[i * dims + 2];

		position_sums[cell_id][0] += px;
		position_sums[cell_id][1] += py;
		position_sums[cell_id][2] += pz;
		position_counts[cell_id]++;

		// Velocity accumulation
		real_t const vx = mech_data.velocities[i * dims];
		real_t const vy = mech_data.velocities[i * dims + 1];
		real_t const vz = mech_data.velocities[i * dims + 2];

		velocity_sums[cell_id][0] += vx;
		velocity_sums[cell_id][1] += vy;
		velocity_sums[cell_id][2] += vz;
		velocity_counts[cell_id]++;

		// Pressure: sum of |force|
		real_t const fx = mech_data.forces[i * dims];
		real_t const fy = mech_data.forces[i * dims + 1];
		real_t const fz = mech_data.forces[i * dims + 2];
		real_t const force_magnitude = std::sqrt(fx * fx + fy * fy + fz * fz);
		e.cells.add_pressure(cell_id, agent_type, force_magnitude);
	}

	// Second pass: compute averages and derived quantities
	for (const auto& [cell_id, sum] : position_sums)
	{
		// Position: average of all agents in the cell
		auto const n = static_cast<real_t>(position_counts[cell_id]);
		e.cells.positions[cell_id] = { sum[0] / n, sum[1] / n, sum[2] / n };
	}

	for (const auto& [cell_id, sum] : velocity_sums)
	{
		auto const n = static_cast<real_t>(velocity_counts[cell_id]);
		real_t const vx = sum[0] / n;
		real_t const vy = sum[1] / n;
		real_t const vz = sum[2] / n;

		e.cells.velocities[cell_id] = { vx, vy, vz };

		// Speed (velocity magnitude)
		real_t const speed = std::sqrt(vx * vx + vy * vy + vz * vz);
		e.cells.speeds[cell_id] = speed;

		// Motility direction (normalized velocity)
		if (speed > 1e-10)
		{
			e.cells.motility_directions[cell_id] = { vx / speed, vy / speed, vz / speed };
		}
		else
		{
			e.cells.motility_directions[cell_id] = { 0.0, 0.0, 0.0 };
		}
	}

	// Third pass: build cell neighbor lists from agent neighbors
	for (index_t i = 0; i < count; ++i)
	{
		index_t const cell_id_i = mech_data.cell_ids[i];
		if (std::cmp_equal(cell_id_i ,-1)))
			continue;

		for (index_t const j : mech_data.neighbors[i])
		{
			index_t const cell_id_j = mech_data.cell_ids[j];
			if (std::cmp_equal(cell_id_j ,-1)))
				continue;

			// Different cells that touch are neighbors
			if (cell_id_i != cell_id_j)
			{
				e.cells.neighbor_cells[cell_id_i].insert(cell_id_j);
			}
		}
	}
}

void openmp_solver::update_motility(environment& e) { m_solver_.update_motility(e); }

void openmp_solver::update_basement_membrane_interactions(environment& e) { bm_solver_.update_interactions(e); }

void openmp_solver::update_spring_attachments(environment& e) { s_solver_.update_spring_attachments(e); }

void openmp_solver::update_positions(environment& e) { p_solver_.update_positions(e); }

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
