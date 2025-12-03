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
	pugi::xml_node child = get_required_child(node, tag_name);
	return static_cast<real_t>(child.text().as_double());
}

bool parse_bool(const pugi::xml_node& node, const char* tag_name)
{
	pugi::xml_node child = get_required_child(node, tag_name);
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

} // namespace

mechanics_config parse_simulation_parameters(const std::filesystem::path& config_file)
{
	if (!std::filesystem::exists(config_file))
	{
		throw std::runtime_error("Configuration file not found: " + config_file.string());
	}

	pugi::xml_document doc;
	if (pugi::xml_parse_result result = doc.load_file(config_file.string().c_str()); !result)
	{
		throw std::runtime_error("Failed to parse XML file: " + config_file.string() + " - " + result.description());
	}

	pugi::xml_node root = doc.child("PhysiCell_settings");
	if (!root)
	{
		throw std::runtime_error("Root <PhysiCell_settings> tag not found in " + config_file.string());
	}

	mechanics_config config {};
	SimulationParameters& params = config.parameters;

	if (pugi::xml_node domain_node = root.child("domain"); domain_node)
	{
		params.is_2D = parse_bool(domain_node, "use_2D");
	}
	else
	{
		params.is_2D = false;
	}

	std::unordered_map<std::string, std::size_t> substrate_index;
	if (pugi::xml_node microenv_node = root.child("microenvironment_setup"); microenv_node)
	{
		for (pugi::xml_node variable_node = microenv_node.child("variable"); variable_node;
			 variable_node = variable_node.next_sibling("variable"))
		{
			std::string name = variable_node.attribute("name").as_string();
			std::size_t idx = config.substrate_names.size();
			config.substrate_names.push_back(name);
			substrate_index.emplace(name, idx);
		}
	}

	std::size_t substrate_count = config.substrate_names.size();

	pugi::xml_node cell_defs_node = get_required_child(root, "cell_definitions");
	std::unordered_map<std::string, std::size_t> cell_name_to_id;
	std::size_t max_id = 0;
	for (pugi::xml_node cell_def = cell_defs_node.child("cell_definition"); cell_def;
		 cell_def = cell_def.next_sibling("cell_definition"))
	{
		if (!cell_def.attribute("ID"))
		{
			throw std::runtime_error("<cell_definition> missing ID attribute");
		}

		std::size_t id = cell_def.attribute("ID").as_uint();
		std::string name = cell_def.attribute("name").as_string();

		if (config.cell_definition_names.size() <= id)
		{
			config.cell_definition_names.resize(id + 1);
		}

		config.cell_definition_names[id] = name;
		cell_name_to_id[name] = id;
		max_id = std::max(max_id, id);
	}

	if (config.cell_definition_names.empty())
	{
		throw std::runtime_error("No <cell_definition> entries found in " + config_file.string());
	}

	std::size_t cell_count = max_id + 1;

	ensure_size(params.is_movable, cell_count, false);
	ensure_size(params.cell_cell_adhesion_strength, cell_count, 0.0);
	ensure_size(params.cell_cell_repulsion_strength, cell_count, 0.0);
	ensure_size(params.relative_maximum_adhesion_distance, cell_count, 0.0);
	ensure_size(params.set_relative_maximum_adhesion_distance, cell_count, 0.0);
	ensure_size(params.set_absolute_maximum_adhesion_distance, cell_count, 0.0);
	ensure_size(params.maximum_number_of_attachments, cell_count, index_t(0));
	ensure_size(params.cell_BM_adhesion_strength, cell_count, 0.0);
	ensure_size(params.cell_BM_repulsion_strength, cell_count, 0.0);
	ensure_size(params.attachment_elastic_coefficient, cell_count, 0.0);
	ensure_size(params.attachment_rate, cell_count, 0.0);
	ensure_size(params.detachment_rate, cell_count, 0.0);
	ensure_size(params.motility_speed, cell_count, 0.0);
	ensure_size(params.motility_persistence_time, cell_count, 0.0);
	ensure_size(params.motility_bias, cell_count, 0.0);
	ensure_size(params.normalize_each_gradient, cell_count, false);

	ensure_size(params.cell_adhesion_affinity, cell_count * cell_count, 0.0);
	ensure_size(params.chemotaxis_sensitivity, cell_count * substrate_count, 0.0);
	ensure_size(params.chemotaxis_enabled, cell_count * substrate_count, false);
	ensure_size(params.chemotaxis_advanced_enabled, cell_count * substrate_count, 0.0);

	params.basic_motility_enabled = false;
	params.advanced_motility_enabled = false;

	for (pugi::xml_node cell_def = cell_defs_node.child("cell_definition"); cell_def;
		 cell_def = cell_def.next_sibling("cell_definition"))
	{
		std::size_t id = cell_def.attribute("ID").as_uint();
		pugi::xml_node phenotype_node = get_required_child(cell_def, "phenotype");

		pugi::xml_node mechanics_node = get_required_child(phenotype_node, "mechanics");
		params.cell_cell_adhesion_strength[id] = parse_real(mechanics_node, "cell_cell_adhesion_strength");
		params.cell_cell_repulsion_strength[id] = parse_real(mechanics_node, "cell_cell_repulsion_strength");
		params.relative_maximum_adhesion_distance[id] =
			parse_real(mechanics_node, "relative_maximum_adhesion_distance");

		if (pugi::xml_node affinities_node = mechanics_node.child("cell_adhesion_affinities"); affinities_node)
		{
			for (pugi::xml_node affinity = affinities_node.child("cell_adhesion_affinity"); affinity;
				 affinity = affinity.next_sibling("cell_adhesion_affinity"))
			{
				std::string other_name = affinity.attribute("name").as_string();
				auto it = cell_name_to_id.find(other_name);
				if (it == cell_name_to_id.end())
				{
					throw std::runtime_error("Unknown cell type in <cell_adhesion_affinity>: " + other_name);
				}

				params.cell_adhesion_affinity[id * cell_count + it->second] =
					static_cast<real_t>(affinity.text().as_double());
			}
		}

		if (pugi::xml_node bm_adh = mechanics_node.child("cell_BM_adhesion_strength"); bm_adh)
		{
			params.cell_BM_adhesion_strength[id] = static_cast<real_t>(bm_adh.text().as_double());
		}

		if (pugi::xml_node bm_rep = mechanics_node.child("cell_BM_repulsion_strength"); bm_rep)
		{
			params.cell_BM_repulsion_strength[id] = static_cast<real_t>(bm_rep.text().as_double());
		}

		if (pugi::xml_node options_node = mechanics_node.child("options"); options_node)
		{
			if (pugi::xml_node rel = options_node.child("set_relative_equilibrium_distance"); rel)
			{
				if (rel.attribute("enabled").as_bool())
				{
					params.set_relative_maximum_adhesion_distance[id] = static_cast<real_t>(rel.text().as_double());
				}
			}

			if (pugi::xml_node abs_node = options_node.child("set_absolute_equilibrium_distance"); abs_node)
			{
				if (abs_node.attribute("enabled").as_bool())
				{
					params.set_absolute_maximum_adhesion_distance[id] =
						static_cast<real_t>(abs_node.text().as_double());
				}
			}
		}

		params.attachment_elastic_coefficient[id] = parse_real(mechanics_node, "attachment_elastic_constant");
		params.attachment_rate[id] = parse_real(mechanics_node, "attachment_rate");
		params.detachment_rate[id] = parse_real(mechanics_node, "detachment_rate");
		if (pugi::xml_node max_attach = mechanics_node.child("maximum_number_of_attachments"); max_attach)
		{
			params.maximum_number_of_attachments[id] = max_attach.text().as_int(0);
		}

		pugi::xml_node motility_node = get_required_child(phenotype_node, "motility");
		params.motility_speed[id] = parse_real(motility_node, "speed");
		params.motility_persistence_time[id] = parse_real(motility_node, "persistence_time");
		params.motility_bias[id] = parse_real(motility_node, "migration_bias");

		if (pugi::xml_node options_node = motility_node.child("options"); options_node)
		{
			bool enabled = options_node.child("enabled") && options_node.child("enabled").text().as_bool();
			bool use_2d = options_node.child("use_2D") && options_node.child("use_2D").text().as_bool();

			params.is_movable[id] = enabled;
			params.basic_motility_enabled = params.basic_motility_enabled || enabled;
			params.is_2D = params.is_2D || use_2d;

			if (pugi::xml_node chemotaxis_node = options_node.child("chemotaxis"); chemotaxis_node)
			{
				bool chem_enabled =
					chemotaxis_node.child("enabled") && chemotaxis_node.child("enabled").text().as_bool();
				std::string substrate = chemotaxis_node.child("substrate").text().as_string();
				auto it = substrate_index.find(substrate);

				if (chem_enabled && it == substrate_index.end())
				{
					throw std::runtime_error("Unknown substrate in <chemotaxis>: " + substrate);
				}

				if (it != substrate_index.end())
				{
					params.chemotaxis_enabled[id * substrate_count + it->second] = chem_enabled;
					if (chem_enabled)
					{
						params.chemotaxis_sensitivity[id * substrate_count + it->second] =
							static_cast<real_t>(chemotaxis_node.child("direction").text().as_double());
					}
				}
			}

			if (pugi::xml_node advanced_node = options_node.child("advanced_chemotaxis"); advanced_node)
			{
				bool advanced_enabled =
					advanced_node.child("enabled") && advanced_node.child("enabled").text().as_bool();
				params.advanced_motility_enabled = params.advanced_motility_enabled || advanced_enabled;

				bool normalize = advanced_node.child("normalize_each_gradient")
								 && advanced_node.child("normalize_each_gradient").text().as_bool();
				params.normalize_each_gradient[id] = advanced_enabled ? normalize : false;

				if (pugi::xml_node sensitivities_node = advanced_node.child("chemotactic_sensitivities");
					sensitivities_node)
				{
					for (pugi::xml_node sens = sensitivities_node.child("chemotactic_sensitivity"); sens;
						 sens = sens.next_sibling("chemotactic_sensitivity"))
					{
						std::string substrate = sens.attribute("substrate").as_string();
						auto it = substrate_index.find(substrate);

						if (advanced_enabled && it == substrate_index.end())
						{
							throw std::runtime_error("Unknown substrate in <advanced_chemotaxis>: " + substrate);
						}

						if (it != substrate_index.end())
						{
							real_t value = static_cast<real_t>(sens.text().as_double());
							if (advanced_enabled)
							{
								params.chemotaxis_advanced_enabled[id * substrate_count + it->second] = value;
								params.chemotaxis_enabled[id * substrate_count + it->second] = true;
								params.chemotaxis_sensitivity[id * substrate_count + it->second] = value;
							}
						}
					}
				}
			}
		}
	}

	return config;
}

} // namespace physicore::mechanics::physicell
