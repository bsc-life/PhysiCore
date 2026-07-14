#include "vtk_serializer.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

#include "microenvironment.h"

using namespace physicore::biofvm;

namespace {
using physicore::real_t;
constexpr const char* vtk_real_type_name()
{
	return std::is_same_v<real_t, float> ? "Float32" : "Float64";
}
} // namespace

vtk_serializer::vtk_serializer(std::string_view output_dir, microenvironment& m)
	: vtk_serializer_base(output_dir, "vtk_microenvironment", "microenvironment.pvd"),
	  substrates_count_(m.substrates_count),
	  substrate_names_(m.substrates_names)
{
	auto x0 = static_cast<int>(m.mesh.bounding_box_mins[0] / static_cast<sindex_t>(m.mesh.voxel_shape[0]));
	auto y0 = static_cast<int>(m.mesh.bounding_box_mins[1] / static_cast<sindex_t>(m.mesh.voxel_shape[1]));
	auto z0 = static_cast<int>(m.mesh.bounding_box_mins[2] / static_cast<sindex_t>(m.mesh.voxel_shape[2]));

	extent_ = { x0,
			   x0 + static_cast<int>(m.mesh.grid_shape[0]),
			   y0,
			   y0 + static_cast<int>(m.mesh.grid_shape[1]),
			   (m.mesh.dims == 3) ? z0 : 0,
			   (m.mesh.dims == 3) ? z0 + static_cast<int>(m.mesh.grid_shape[2]) : 0 };

	spacing_ = { static_cast<double>(m.mesh.voxel_shape[0]),
				 static_cast<double>(m.mesh.voxel_shape[1]),
				 (m.mesh.dims == 3) ? static_cast<double>(m.mesh.voxel_shape[2]) : 0.0 };
}

void vtk_serializer::serialize(const microenvironment& m, real_t current_time)
{
	std::ostringstream name_ss;
	name_ss << "microenvironment_" << std::setw(6) << std::setfill('0') << iteration << ".vti";
	const auto file_name = name_ss.str();
	const auto file_path = std::filesystem::path(vtks_dir) / file_name;

	const auto [e0, e1, e2, e3, e4, e5] = extent_;
	const auto [sx, sy, sz] = spacing_;
	const auto n_voxels = m.mesh.voxel_count();

	std::ofstream out(file_path);
	out << std::setprecision(std::numeric_limits<real_t>::max_digits10);

	out << "<?xml version=\"1.0\"?>\n"
		<< "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
		<< "  <ImageData WholeExtent=\"" << e0 << " " << e1 << " " << e2 << " " << e3 << " " << e4 << " " << e5
		<< "\" Origin=\"0 0 0\" Spacing=\"" << sx << " " << sy << " " << sz << "\">\n"
		<< "    <Piece Extent=\"" << e0 << " " << e1 << " " << e2 << " " << e3 << " " << e4 << " " << e5 << "\">\n"
		<< "      <CellData>\n";

	for (index_t s = 0; s < substrates_count_; ++s)
	{
		out << "        <DataArray type=\"" << vtk_real_type_name() << "\" Name=\"" << substrate_names_[s]
			<< "\" NumberOfComponents=\"1\" NumberOfTuples=\"" << n_voxels << "\" format=\"ascii\">\n          ";
		for (index_t z = 0; z < m.mesh.grid_shape[2]; ++z)
			for (index_t y = 0; y < m.mesh.grid_shape[1]; ++y)
				for (index_t x = 0; x < m.mesh.grid_shape[0]; ++x)
					out << m.get_substrate_density(s, x, y, z) << " ";
		out << "\n        </DataArray>\n";
	}

	out << "      </CellData>\n    </Piece>\n  </ImageData>\n</VTKFile>\n";

	append_to_pvd(file_name, current_time);
	iteration++;
}
