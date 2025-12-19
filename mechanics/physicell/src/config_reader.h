#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "mechanical_parameters.h"

namespace physicore::mechanics::physicell {

struct domain_config
{
	real_t x_min;
	real_t x_max;
	real_t y_min;
	real_t y_max;
	real_t z_min;
	real_t z_max;
	real_t dx;
	real_t dy;
	real_t dz;
	bool use_2D;
};

struct overall_config
{
	real_t max_time;
	std::string time_units;
	std::string space_units;
	real_t dt_mechanics;
};

struct mechanics_config
{
	domain_config domain;
	overall_config overall;
	std::vector<mechanical_parameters> cell_types; // One entry per cell type, indexed by ID
	bool is_2D;									   // Global domain setting
};

/**
 * @brief Parse mechanics-related settings from a PhysiCell_settings.xml file.
 *
 * Reads the mechanics / motility sections for each <cell_definition> and creates
 * mechanical_parameters entries indexed by the ID declared in the XML.
 *
 * @param config_file Path to PhysiCell_settings.xml
 * @return mechanics_config containing cell type parameters and substrate names
 * @throws std::runtime_error if the file is missing, malformed, or references unknown cell types / substrates
 */
mechanics_config parse_simulation_parameters(const std::filesystem::path& config_file);

} // namespace physicore::mechanics::physicell
