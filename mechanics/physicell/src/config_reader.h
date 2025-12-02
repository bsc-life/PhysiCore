#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "simulations_parameters.h"

namespace physicore::mechanics::physicell {

struct mechanics_config
{
	SimulationParameters<> parameters;
	std::vector<std::string> cell_definition_names; // Indexed by cell_definition ID
	std::vector<std::string> substrate_names;		   // Order taken from <microenvironment_setup>
};

/**
 * @brief Parse mechanics-related settings from a PhysiCell_settings.xml file.
 *
 * Reads the mechanics / motility sections for each <cell_definition> and fills
 * SimulationParameters using the IDs declared in the XML as indices.
 *
 * @param config_file Path to PhysiCell_settings.xml
 * @return mechanics_config containing SimulationParameters and name lookups
 * @throws std::runtime_error if the file is missing, malformed, or references unknown cell types / substrates
 */
mechanics_config parse_simulation_parameters(const std::filesystem::path& config_file);

} // namespace physicore::mechanics::physicell
