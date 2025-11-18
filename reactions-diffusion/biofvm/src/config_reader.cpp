#include "config_reader.h"

#include <pugixml.hpp>
#include <stdexcept>
#include <string>

namespace physicore {

namespace {

// Helper function to get required child node
pugi::xml_node get_required_child(const pugi::xml_node& parent, const char* name)
{
	pugi::xml_node child = parent.child(name);
	if (!child)
	{
		throw std::runtime_error(std::string("Required XML tag <") + name + "> not found under <" + parent.name()
								 + ">");
	}
	return child;
}

// Helper function to parse text content as real_t
real_t parse_real(const pugi::xml_node& node, const char* tag_name)
{
	pugi::xml_node child = get_required_child(node, tag_name);
	return static_cast<real_t>(child.text().as_double());
}

// Helper function to parse text content as bool
bool parse_bool(const pugi::xml_node& node, const char* tag_name)
{
	pugi::xml_node child = get_required_child(node, tag_name);
	return child.text().as_bool();
}

// Helper function to parse text content as string
std::string parse_string(const pugi::xml_node& node, const char* tag_name)
{
	pugi::xml_node child = get_required_child(node, tag_name);
	return child.text().as_string();
}

// Parse <domain> tag
domain_config parse_domain(const pugi::xml_node& domain_node)
{
	domain_config config;
	config.x_min = parse_real(domain_node, "x_min");
	config.x_max = parse_real(domain_node, "x_max");
	config.y_min = parse_real(domain_node, "y_min");
	config.y_max = parse_real(domain_node, "y_max");
	config.z_min = parse_real(domain_node, "z_min");
	config.z_max = parse_real(domain_node, "z_max");
	config.dx = parse_real(domain_node, "dx");
	config.dy = parse_real(domain_node, "dy");
	config.dz = parse_real(domain_node, "dz");
	config.use_2D = parse_bool(domain_node, "use_2D");
	return config;
}

// Parse <overall> tag
overall_config parse_overall(const pugi::xml_node& overall_node)
{
	overall_config config;
	config.max_time = parse_real(overall_node, "max_time");
	config.time_units = parse_string(overall_node, "time_units");
	config.space_units = parse_string(overall_node, "space_units");
	config.dt_diffusion = parse_real(overall_node, "dt_diffusion");
	config.dt_mechanics = parse_real(overall_node, "dt_mechanics");
	config.dt_phenotype = parse_real(overall_node, "dt_phenotype");
	return config;
}

// Parse Dirichlet boundary conditions for a variable
dirichlet_boundary_config parse_dirichlet_options(const pugi::xml_node& variable_node)
{
	dirichlet_boundary_config config;

	// Initialize with defaults (disabled, zero values)
	config.mins_values.fill(0.0);
	config.maxs_values.fill(0.0);
	config.mins_conditions.fill(false);
	config.maxs_conditions.fill(false);

	{
		// Check for legacy Dirichlet_boundary_condition
		pugi::xml_node boundary_node = variable_node.child("Dirichlet_boundary_condition");
		if (boundary_node)
		{
			bool enabled = boundary_node.attribute("enabled").as_bool();
			real_t value = static_cast<real_t>(boundary_node.text().as_double());
			if (enabled)
			{
				// Apply to all boundaries
				config.mins_values.fill(value);
				config.maxs_values.fill(value);
				config.mins_conditions.fill(true);
				config.maxs_conditions.fill(true);
			}
		}
	}

	pugi::xml_node options_node = variable_node.child("Dirichlet_options");
	if (!options_node)
	{
		return config;
	}

	// Parse individual boundary values
	const char* boundary_ids[] = { "xmin", "xmax", "ymin", "ymax", "zmin", "zmax" };
	for (int i = 0; i < 6; ++i)
	{
		pugi::xml_node boundary_value = options_node.find_child_by_attribute("boundary_value", "ID", boundary_ids[i]);
		if (boundary_value)
		{
			real_t value = static_cast<real_t>(boundary_value.text().as_double());
			bool enabled = boundary_value.attribute("enabled").as_bool();

			if (i % 2 == 0)
			{
				// mins: xmin=0, ymin=2, zmin=4
				int axis = i / 2;
				config.mins_values[axis] = value;
				config.mins_conditions[axis] = enabled;
			}
			else
			{
				// maxs: xmax=1, ymax=3, zmax=5
				int axis = i / 2;
				config.maxs_values[axis] = value;
				config.maxs_conditions[axis] = enabled;
			}
		}
	}

	return config;
}

// Parse a single <variable> tag
variable_config parse_variable(const pugi::xml_node& variable_node)
{
	variable_config config;

	// Parse attributes
	config.name = variable_node.attribute("name").as_string();
	config.units = variable_node.attribute("units").as_string();
	config.id = variable_node.attribute("ID").as_uint();

	// Parse physical parameters
	pugi::xml_node param_set = get_required_child(variable_node, "physical_parameter_set");
	config.diffusion_coefficient = parse_real(param_set, "diffusion_coefficient");
	config.decay_rate = parse_real(param_set, "decay_rate");

	// Parse initial condition
	config.initial_condition = parse_real(variable_node, "initial_condition");

	// Parse Dirichlet boundary conditions
	config.boundary_conditions = parse_dirichlet_options(variable_node);

	return config;
}

// Parse <microenvironment_setup> tag
microenvironment_config parse_microenvironment(const pugi::xml_node& microenv_node)
{
	microenvironment_config config;

	// Parse all variables
	for (pugi::xml_node variable_node = microenv_node.child("variable"); variable_node;
		 variable_node = variable_node.next_sibling("variable"))
	{
		config.variables.push_back(parse_variable(variable_node));
	}

	if (config.variables.empty())
	{
		throw std::runtime_error("No <variable> tags found in <microenvironment_setup>");
	}

	// Parse options
	pugi::xml_node options_node = microenv_node.child("options");
	if (options_node)
	{
		pugi::xml_node grad_node = options_node.child("calculate_gradients");
		config.calculate_gradients = grad_node ? grad_node.text().as_bool() : false;

		pugi::xml_node track_node = options_node.child("track_internalized_substrates_in_each_agent");
		config.track_internalized_substrates = track_node ? track_node.text().as_bool() : false;
	}
	else
	{
		config.calculate_gradients = false;
		config.track_internalized_substrates = false;
	}

	return config;
}

} // namespace

physicell_config parse_physicell_config(const std::filesystem::path& config_file)
{
	// Check file exists
	if (!std::filesystem::exists(config_file))
	{
		throw std::runtime_error("Configuration file not found: " + config_file.string());
	}

	// Load XML document
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file(config_file.string().c_str());

	if (!result)
	{
		throw std::runtime_error("Failed to parse XML file: " + config_file.string() + " - " + result.description());
	}

	// Get root node
	pugi::xml_node root = doc.child("PhysiCell_settings");
	if (!root)
	{
		throw std::runtime_error("Root <PhysiCell_settings> tag not found in " + config_file.string());
	}

	physicell_config config;

	// Parse required sections
	pugi::xml_node domain_node = get_required_child(root, "domain");
	config.domain = parse_domain(domain_node);

	pugi::xml_node overall_node = get_required_child(root, "overall");
	config.overall = parse_overall(overall_node);

	pugi::xml_node microenv_node = get_required_child(root, "microenvironment_setup");
	config.microenvironment = parse_microenvironment(microenv_node);

	return config;
}

} // namespace physicore
