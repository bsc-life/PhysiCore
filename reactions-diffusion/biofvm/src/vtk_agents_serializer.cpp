#include "vtk_agents_serializer.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

#include "agent_container.h"
#include "microenvironment.h"

using namespace physicore::biofvm;

namespace {
using physicore::real_t;
constexpr const char* vtk_real_type_name()
{
	return std::is_same_v<real_t, float> ? "Float32" : "Float64";
}
} // namespace

vtk_agents_serializer::vtk_agents_serializer(std::string_view output_dir, const microenvironment& m)
	: vtk_serializer_base(output_dir, "vtk_agents", "agents.pvd"),
	  substrate_count_(m.substrates_count),
	  substrate_names_(m.substrates_names)
{
}

void vtk_agents_serializer::serialize(const microenvironment& m, real_t current_time)
{
	const auto& biofvm_data = retrieve_agent_data(*m.agents);
	const auto& base_data = biofvm_data.base_data;
	const index_t agent_count = biofvm_data.agents_count;

	std::ostringstream name_ss;
	name_ss << "agents_" << std::setw(6) << std::setfill('0') << iteration << ".vtu";
	const auto file_name = name_ss.str();
	const auto file_path = std::filesystem::path(vtks_dir) / file_name;

	std::ofstream out(file_path);
	out << std::setprecision(std::numeric_limits<real_t>::max_digits10);

	out << "<?xml version=\"1.0\"?>\n"
		<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
		<< "  <UnstructuredGrid>\n"
		<< "    <Piece NumberOfPoints=\"" << agent_count << "\" NumberOfCells=\"" << agent_count << "\">\n";

	// Points
	out << "      <Points>\n"
		<< "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n"
		<< "          ";
	for (index_t i = 0; i < agent_count; ++i)
	{
		std::array<double, 3> pos = { 0.0, 0.0, 0.0 };
		for (index_t d = 0; d < base_data.dims && d < 3; ++d)
			pos[d] = base_data.positions[i * base_data.dims + d];
		out << pos[0] << " " << pos[1] << " " << pos[2] << " ";
	}
	out << "\n        </DataArray>\n      </Points>\n";

	// Cells — one VTK_VERTEX (type=1) per agent
	out << "      <Cells>\n"
		<< "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n          ";
	for (index_t i = 0; i < agent_count; ++i) out << i << " ";
	out << "\n        </DataArray>\n"
		<< "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n          ";
	for (index_t i = 0; i < agent_count; ++i) out << (i + 1) << " ";
	out << "\n        </DataArray>\n"
		<< "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ";
	for (index_t i = 0; i < agent_count; ++i) out << "1 ";
	out << "\n        </DataArray>\n      </Cells>\n";

	// PointData arrays
	out << "      <PointData>\n";

	auto write_array = [&](const std::string& name, auto get_value) {
		out << "        <DataArray type=\"" << vtk_real_type_name() << "\" Name=\"" << name
			<< "\" NumberOfComponents=\"1\" NumberOfTuples=\"" << agent_count << "\" format=\"ascii\">\n          ";
		for (index_t i = 0; i < agent_count; ++i) out << get_value(i) << " ";
		out << "\n        </DataArray>\n";
	};

	write_array("volume", [&](index_t i) { return biofvm_data.volumes[i]; });

	for (index_t s = 0; s < substrate_count_; ++s)
	{
		const auto& sname = substrate_names_[s];
		write_array(sname + "_secretion_rate",
					[&](index_t i) { return biofvm_data.secretion_rates[i * substrate_count_ + s]; });
		write_array(sname + "_saturation_density",
					[&](index_t i) { return biofvm_data.saturation_densities[i * substrate_count_ + s]; });
		write_array(sname + "_uptake_rate",
					[&](index_t i) { return biofvm_data.uptake_rates[i * substrate_count_ + s]; });
		write_array(sname + "_net_export_rate",
					[&](index_t i) { return biofvm_data.net_export_rates[i * substrate_count_ + s]; });
		write_array(sname + "_internalized_substrate",
					[&](index_t i) { return biofvm_data.internalized_substrates[i * substrate_count_ + s]; });
		write_array(sname + "_fraction_released_at_death",
					[&](index_t i) { return biofvm_data.fraction_released_at_death[i * substrate_count_ + s]; });
		write_array(sname + "_fraction_transferred_when_ingested",
					[&](index_t i) { return biofvm_data.fraction_transferred_when_ingested[i * substrate_count_ + s]; });
	}

	out << "      </PointData>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";

	append_to_pvd(file_name, current_time);
	iteration++;
}
