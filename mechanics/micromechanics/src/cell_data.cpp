#include "micromechanics/cell_data.h"

namespace physicore::mechanics::micromechanics {

void cell_data::clear()
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

real_t cell_data::get_pressure(index_t cell_id, std::uint8_t agent_type) const
{
	auto it = compartment_pressures.find({ cell_id, agent_type });
	return (it != compartment_pressures.end()) ? it->second : 0.0;
}

real_t cell_data::get_total_cell_pressure(index_t cell_id) const
{
	real_t total = 0.0;
	for (const auto& [key, pressure] : compartment_pressures)
	{
		if (key.first == cell_id)
			total += pressure;
	}
	return total;
}

void cell_data::add_pressure(index_t cell_id, std::uint8_t agent_type, real_t delta)
{
	compartment_pressures[{ cell_id, agent_type }] += delta;
}

index_t cell_data::get_compartment_count(index_t cell_id, std::uint8_t agent_type) const
{
	auto it = compartment_counts.find({ cell_id, agent_type });
	return (it != compartment_counts.end()) ? it->second : 0;
}

index_t cell_data::get_total_agent_count(index_t cell_id) const
{
	index_t total = 0;
	for (const auto& [key, count] : compartment_counts)
	{
		if (key.first == cell_id)
			total += count;
	}
	return total;
}

} // namespace physicore::mechanics::micromechanics
