#pragma once

#include <array>
#include <filesystem>
#include <string>
#include <vector>

#include "types.h"

namespace physicore {

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
	real_t dt_diffusion;
	real_t dt_mechanics; // Stored for future use by mechanics module
	real_t dt_phenotype; // Stored for future use by phenotype module
};

struct dirichlet_boundary_config
{
	std::array<real_t, 3> mins_values;	 // [xmin, ymin, zmin]
	std::array<real_t, 3> maxs_values;	 // [xmax, ymax, zmax]
	std::array<bool, 3> mins_conditions; // [xmin_enabled, ymin_enabled, zmin_enabled]
	std::array<bool, 3> maxs_conditions; // [xmax_enabled, ymax_enabled, zmax_enabled]
};

struct variable_config
{
	std::string name;
	std::string units;
	index_t id;
	real_t diffusion_coefficient;
	real_t decay_rate;
	real_t initial_condition;
	dirichlet_boundary_config boundary_conditions;
};

struct microenvironment_config
{
	std::vector<variable_config> variables;
	bool calculate_gradients;
	bool track_internalized_substrates;
};

struct physicell_config
{
	domain_config domain;
	overall_config overall;
	microenvironment_config microenvironment;
};

/**
 * @brief Parse PhysiCell_settings.xml file and extract configuration.
 *
 * This function performs structural validation:
 * - Verifies file exists and is readable
 * - Checks required XML tags are present (<domain>, <overall>, <microenvironment_setup>)
 * - Validates XML structure and parsability
 *
 * Domain validation (e.g., positive values, consistent bounds) is deferred to
 * the builder or constructor that consumes this configuration.
 *
 * @param config_file Path to PhysiCell_settings.xml file
 * @return Parsed configuration structures
 * @throws std::runtime_error if file cannot be read or XML is malformed/incomplete
 */
physicell_config parse_physicell_config(const std::filesystem::path& config_file);

} // namespace physicore
