#include "config_reader.h"

#include <algorithm>
#include <pugixml.hpp>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace physicore::mechanics::physicell {

namespace {

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

real_t parse_real(const pugi::xml_node& node, const char* tag_name)
{
	const pugi::xml_node child = get_required_child(node, tag_name);
	return static_cast<real_t>(child.text().as_double());
}

bool parse_bool(const pugi::xml_node& node, const char* tag_name)
{
	const pugi::xml_node child = get_required_child(node, tag_name);
	return child.text().as_bool();
}

template <typename Value>
void ensure_size(std::vector<Value>& vec, std::size_t size, const Value& init)
{
	vec.assign(size, init);
}

template <typename Value>
void ensure_matrix(std::vector<std::vector<Value>>& matrix, std::size_t rows, std::size_t cols, const Value& init)
{
	matrix.assign(rows, std::vector<Value>(cols, init));
}

// Parse cell adhesion affinities for a cell type
void parse_cell_affinities(const pugi::xml_node& mechanics_node, mechanical_parameters& params,
						   const std::unordered_map<std::string, std::size_t>& cell_name_to_id)
{
	if (const pugi::xml_node affinities_node = mechanics_node.child("cell_adhesion_affinities"); affinities_node)
	{
		for (pugi::xml_node affinity = affinities_node.child("cell_adhesion_affinity"); affinity;
			 affinity = affinity.next_sibling("cell_adhesion_affinity"))
		{
			const std::string other_name = affinity.attribute("name").as_string();
			auto it = cell_name_to_id.find(other_name);
			if (it == cell_name_to_id.end())
			{
				throw std::runtime_error("Unknown cell type in <cell_adhesion_affinity>: " + other_name);
			}
			params.cell_adhesion_affinity[it->second] = static_cast<real_t>(affinity.text().as_double());
		}
	}
}

// Parse mechanics parameters from XML
void parse_mechanics_params(const pugi::xml_node& mechanics_node, mechanical_parameters& params,
							const std::unordered_map<std::string, std::size_t>& cell_name_to_id)
{
	params.cell_cell_adhesion_strength = parse_real(mechanics_node, "cell_cell_adhesion_strength");
	params.cell_cell_repulsion_strength = parse_real(mechanics_node, "cell_cell_repulsion_strength");
	params.relative_maximum_adhesion_distance = parse_real(mechanics_node, "relative_maximum_adhesion_distance");

	parse_cell_affinities(mechanics_node, params, cell_name_to_id);

	// Optional BM adhesion/repulsion
	if (const pugi::xml_node bm_adh = mechanics_node.child("cell_BM_adhesion_strength"); bm_adh)
	{
		params.cell_BM_adhesion_strength = static_cast<real_t>(bm_adh.text().as_double());
	}
	if (const pugi::xml_node bm_rep = mechanics_node.child("cell_BM_repulsion_strength"); bm_rep)
	{
		params.cell_BM_repulsion_strength = static_cast<real_t>(bm_rep.text().as_double());
	}

	// Parse equilibrium distance options
	if (const pugi::xml_node options_node = mechanics_node.child("options"); options_node)
	{
		if (const pugi::xml_node rel = options_node.child("set_relative_equilibrium_distance"); rel)
		{
			if (rel.attribute("enabled").as_bool())
			{
				params.set_relative_maximum_adhesion_distance = static_cast<real_t>(rel.text().as_double());
			}
		}
		if (const pugi::xml_node abs_node = options_node.child("set_absolute_equilibrium_distance"); abs_node)
		{
			if (abs_node.attribute("enabled").as_bool())
			{
				params.set_absolute_maximum_adhesion_distance = static_cast<real_t>(abs_node.text().as_double());
			}
		}
	}

	params.attachment_elastic_coefficient = parse_real(mechanics_node, "attachment_elastic_constant");
	params.attachment_rate = parse_real(mechanics_node, "attachment_rate");
	params.detachment_rate = parse_real(mechanics_node, "detachment_rate");
	if (const pugi::xml_node max_attach = mechanics_node.child("maximum_number_of_attachments"); max_attach)
	{
		params.maximum_number_of_attachments = max_attach.text().as_int(0);
	}
}

// Parse chemotaxis options
void parse_chemotaxis_options(const pugi::xml_node& options_node, mechanical_parameters& params,
							  const std::unordered_map<std::string, std::size_t>& substrate_index)
{
	if (const pugi::xml_node chemotaxis_node = options_node.child("chemotaxis"); chemotaxis_node)
	{
		const bool chem_enabled = chemotaxis_node.child("enabled") && chemotaxis_node.child("enabled").text().as_bool();
		const std::string substrate = chemotaxis_node.child("substrate").text().as_string();
		auto it = substrate_index.find(substrate);

		if (chem_enabled && it == substrate_index.end())
		{
			throw std::runtime_error("Unknown substrate in <chemotaxis>: " + substrate);
		}

		if (it != substrate_index.end())
		{
			params.chemotaxis_enabled[it->second] = chem_enabled;
			if (chem_enabled)
			{
				params.chemotaxis_sensitivity[it->second] =
					static_cast<real_t>(chemotaxis_node.child("direction").text().as_double());
			}
		}
	}
}

void parse_advanced_chemotaxis_sensitivities(const pugi::xml_node& sensitivities_node, const bool advanced_enabled,
											 const std::unordered_map<std::string, std::size_t>& substrate_index,
											 mechanical_parameters& params)
{
	for (pugi::xml_node sens = sensitivities_node.child("chemotactic_sensitivity"); sens;
		 sens = sens.next_sibling("chemotactic_sensitivity"))
	{
		const std::string substrate = sens.attribute("substrate").as_string();
		auto it = substrate_index.find(substrate);

		if (advanced_enabled && it == substrate_index.end())
		{
			throw std::runtime_error("Unknown substrate in <advanced_chemotaxis>: " + substrate);
		}
		if (!advanced_enabled || it == substrate_index.end())
		{
			continue;
		}

		const real_t value = static_cast<real_t>(sens.text().as_double());
		params.chemotaxis_advanced_enabled[it->second] = value;
		params.chemotaxis_enabled[it->second] = true;
		params.chemotaxis_sensitivity[it->second] = value;
	}
}

// Parse advanced chemotaxis options
void parse_advanced_chemotaxis_options(const pugi::xml_node& options_node, mechanical_parameters& params,
									   const std::unordered_map<std::string, std::size_t>& substrate_index)
{
	if (const pugi::xml_node advanced_node = options_node.child("advanced_chemotaxis"); advanced_node)
	{
		const bool advanced_enabled = advanced_node.child("enabled") && advanced_node.child("enabled").text().as_bool();
		params.normalize_each_gradient = advanced_enabled
										 && (advanced_node.child("normalize_each_gradient")
											 && advanced_node.child("normalize_each_gradient").text().as_bool());

		if (const pugi::xml_node sensitivities_node = advanced_node.child("chemotactic_sensitivities");
			sensitivities_node)
		{
			parse_advanced_chemotaxis_sensitivities(sensitivities_node, advanced_enabled, substrate_index, params);
		}
	}
}

// Parse motility parameters and chemotaxis settings
void parse_motility_params(const pugi::xml_node& motility_node, mechanical_parameters& params,
						   const std::unordered_map<std::string, std::size_t>& substrate_index)
{
	params.motility_speed = parse_real(motility_node, "speed");
	params.motility_persistence_time = parse_real(motility_node, "persistence_time");
	params.motility_bias = parse_real(motility_node, "migration_bias");

	if (const pugi::xml_node options_node = motility_node.child("options"); options_node)
	{
		const bool enabled = options_node.child("enabled") && options_node.child("enabled").text().as_bool();
		params.is_movable = enabled;

		parse_chemotaxis_options(options_node, params, substrate_index);
		parse_advanced_chemotaxis_options(options_node, params, substrate_index);
	}
}

// Parse a single cell definition
void parse_cell_definition(pugi::xml_node cell_def, mechanical_parameters& params,
						   const std::unordered_map<std::string, std::size_t>& cell_name_to_id,
						   const std::unordered_map<std::string, std::size_t>& substrate_index,
						   std::size_t agent_type_count)
{
	const std::size_t id = cell_def.attribute("ID").as_uint();
	const std::string name = cell_def.attribute("name").as_string();

	params.id = id;
	params.name = name;
	params.cell_adhesion_affinity.assign(agent_type_count, 0.0);
	params.chemotaxis_sensitivity.assign(substrate_index.size(), 0.0);
	params.chemotaxis_enabled.assign(substrate_index.size(), false);
	params.chemotaxis_advanced_enabled.assign(substrate_index.size(), 0.0);

	const pugi::xml_node phenotype_node = get_required_child(cell_def, "phenotype");

	const pugi::xml_node mechanics_node = get_required_child(phenotype_node, "mechanics");
	parse_mechanics_params(mechanics_node, params, cell_name_to_id);

	const pugi::xml_node motility_node = get_required_child(phenotype_node, "motility");
	parse_motility_params(motility_node, params, substrate_index);
}

// Parse <domain> configuration
domain_config parse_domain(const pugi::xml_node& domain_node)
{
	domain_config config {};
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

// Parse <overall> configuration
overall_config parse_overall(const pugi::xml_node& overall_node)
{
	overall_config config;
	config.max_time = parse_real(overall_node, "max_time");
	config.time_units = get_required_child(overall_node, "time_units").text().as_string();
	config.space_units = get_required_child(overall_node, "space_units").text().as_string();
	config.dt_mechanics = parse_real(overall_node, "dt_mechanics");
	return config;
}

} // namespace

mechanics_config parse_simulation_parameters(const std::filesystem::path& config_file)
{
	if (!std::filesystem::exists(config_file))
	{
		throw std::runtime_error("Configuration file not found: " + config_file.string());
	}

	pugi::xml_document doc;
	if (const pugi::xml_parse_result result = doc.load_file(config_file.string().c_str()); !result)
	{
		throw std::runtime_error("Failed to parse XML file: " + config_file.string() + " - " + result.description());
	}

	const pugi::xml_node root = doc.child("PhysiCell_settings");
	if (!root)
	{
		throw std::runtime_error("Root <PhysiCell_settings> tag not found in " + config_file.string());
	}

	mechanics_config config {};
	config.is_2D = false;

	// Parse domain configuration
	if (const pugi::xml_node domain_node = root.child("domain"); domain_node)
	{
		config.domain = parse_domain(domain_node);
		config.is_2D = config.domain.use_2D;
	}

	// Parse overall configuration
	if (const pugi::xml_node overall_node = root.child("overall"); overall_node)
	{
		config.overall = parse_overall(overall_node);
	}

	std::unordered_map<std::string, std::size_t> substrate_index;
	std::size_t idx = 0;
	if (const pugi::xml_node microenv_node = root.child("microenvironment_setup"); microenv_node)
	{
		for (pugi::xml_node variable_node = microenv_node.child("variable"); variable_node;
			 variable_node = variable_node.next_sibling("variable"))
		{
			const std::string name = variable_node.attribute("name").as_string();
			substrate_index.emplace(name, idx);
			idx++;
		}
	}

	const pugi::xml_node cell_defs_node = get_required_child(root, "cell_definitions");
	std::unordered_map<std::string, std::size_t> cell_name_to_id;
	std::size_t next_id = 0;

	// First pass: map cell names to IDs and determine cell count
	for (pugi::xml_node cell_def = cell_defs_node.child("cell_definition"); cell_def;
		 cell_def = cell_def.next_sibling("cell_definition"))
	{
		if (!cell_def.attribute("ID"))
		{
			throw std::runtime_error("<cell_definition> missing ID attribute");
		}

		const std::size_t id = cell_def.attribute("ID").as_uint();
		if (next_id != id)
		{
			throw std::runtime_error("Cell definition IDs must be sequential starting from 0");
		}
		const std::string name = cell_def.attribute("name").as_string();

		cell_name_to_id[name] = id;
		next_id++;
	}

	const std::size_t cell_type_count = next_id;
	config.cell_types.resize(cell_type_count);

	// Second pass: parse each cell definition
	for (pugi::xml_node cell_def = cell_defs_node.child("cell_definition"); cell_def;
		 cell_def = cell_def.next_sibling("cell_definition"))
	{
		const std::size_t id = cell_def.attribute("ID").as_uint();
		parse_cell_definition(cell_def, config.cell_types[id], cell_name_to_id, substrate_index, cell_type_count);
	}

	return config;
}

} // namespace physicore::mechanics::physicell
