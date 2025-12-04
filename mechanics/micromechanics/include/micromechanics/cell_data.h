#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include <common/types.h>

namespace physicore::mechanics::micromechanics {

/**
 * @brief Hash function for (cell_id, agent_type) pairs.
 */
struct cell_compartment_hash
{
	std::size_t operator()(const std::pair<index_t, std::uint8_t>& p) const noexcept
	{
		// Combine cell_id and agent_type into a single hash
		return std::hash<index_t> {}(p.first) ^ (std::hash<std::uint8_t> {}(p.second) << 1);
	}
};

/**
 * @brief Cell instance data - runtime state for each cell.
 *
 * All properties are indexed by cell_id. Properties are computed
 * by aggregating data from agents belonging to each cell.
 */
struct cell_data
{
	// ========== Geometry ==========

	/// Cell positions indexed by cell_id.
	/// Calculated as the average position of all agents in the cell.
	std::unordered_map<index_t, std::array<real_t, 3>> positions;

	/// Cell volumes indexed by cell_id.
	/// Sum of sub-agent volumes (4/3 * pi * r^3 for each agent).
	std::unordered_map<index_t, real_t> volumes;

	// ========== Kinematics ==========

	/// Cell velocity vectors indexed by cell_id.
	/// Average velocity of all agents in the cell.
	std::unordered_map<index_t, std::array<real_t, 3>> velocities;

	/// Cell speed (velocity magnitude) indexed by cell_id.
	std::unordered_map<index_t, real_t> speeds;

	/// Cell motility direction (normalized) indexed by cell_id.
	/// Direction the cell is actively moving toward.
	std::unordered_map<index_t, std::array<real_t, 3>> motility_directions;

	// ========== Mechanics ==========

	/// Compartment pressures indexed by (cell_id, agent_type).
	/// Pressure = sum(|forces|) for agents in each compartment.
	std::unordered_map<std::pair<index_t, std::uint8_t>, real_t, cell_compartment_hash> compartment_pressures;

	// ========== Topology ==========

	/// Cell neighbor lists indexed by cell_id.
	/// Derived from agent neighbors: if agent_i touches agent_j,
	/// then cell[agent_i] and cell[agent_j] are neighbors.
	std::unordered_map<index_t, std::set<index_t>> neighbor_cells;

	/// Compartment agent counts indexed by (cell_id, agent_type).
	/// Number of agents per compartment type.
	std::unordered_map<std::pair<index_t, std::uint8_t>, index_t, cell_compartment_hash> compartment_counts;

	// ========== Methods ==========

	/// Clear all cell data.
	void clear()
	{
		positions.clear();
		volumes.clear();
		velocities.clear();
		speeds.clear();
		motility_directions.clear();
		compartment_pressures.clear();
		neighbor_cells.clear();
		compartment_counts.clear();
	}

	/// Get pressure for a specific cell and compartment.
	real_t get_pressure(index_t cell_id, std::uint8_t agent_type) const
	{
		auto it = compartment_pressures.find({ cell_id, agent_type });
		return (it != compartment_pressures.end()) ? it->second : 0.0;
	}

	/// Get total pressure for a cell (sum across all compartments).
	real_t get_total_cell_pressure(index_t cell_id) const
	{
		real_t total = 0.0;
		for (const auto& [key, pressure] : compartment_pressures)
		{
			if (key.first == cell_id)
				total += pressure;
		}
		return total;
	}

	/// Add to pressure for a specific cell and compartment.
	void add_pressure(index_t cell_id, std::uint8_t agent_type, real_t delta)
	{
		compartment_pressures[{ cell_id, agent_type }] += delta;
	}

	/// Get agent count for a specific cell and compartment.
	index_t get_compartment_count(index_t cell_id, std::uint8_t agent_type) const
	{
		auto it = compartment_counts.find({ cell_id, agent_type });
		return (it != compartment_counts.end()) ? it->second : 0;
	}

	/// Get total agent count for a cell (sum across all compartments).
	index_t get_total_agent_count(index_t cell_id) const
	{
		index_t total = 0;
		for (const auto& [key, count] : compartment_counts)
		{
			if (key.first == cell_id)
				total += count;
		}
		return total;
	}
};

} // namespace physicore::mechanics::micromechanics
