#include <micromechanics/simulation_parameters.h>

namespace physicore::mechanics::micromechanics {

void simulation_parameters::add_interaction(std::uint8_t type_a, std::uint8_t type_b, const interaction_config& config)
{
	interactions[{ type_a, type_b }] = config;
	if (type_a != type_b)
	{
		interactions[{ type_b, type_a }] = config;
	}
}

const interaction_config& simulation_parameters::get_interaction(std::uint8_t type_a, std::uint8_t type_b) const
{
	auto it = interactions.find({ type_a, type_b });
	if (it != interactions.end())
	{
		return it->second;
	}
	return default_interaction;
}

void simulation_parameters::set_single_type_interaction(const interaction_config& config)
{
	interactions[{ 0, 0 }] = config;
}

} // namespace physicore::mechanics::micromechanics
