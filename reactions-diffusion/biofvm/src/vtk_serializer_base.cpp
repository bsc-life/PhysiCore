#include "vtk_serializer_base.h"

#include <filesystem>
#include <fstream>

using namespace physicore::biofvm;

vtk_serializer_base::vtk_serializer_base(std::string_view output_dir, std::string_view vtks_dir_name)
	: iteration(0), output_dir(output_dir), vtks_dir(std::filesystem::path(output_dir) / vtks_dir_name)
{
	std::filesystem::create_directories(vtks_dir);

	pvd_contents = R"(<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>
)";
}

void vtk_serializer_base::append_to_pvd(std::string_view vtk_file_name)
{
	auto file_path = std::filesystem::path(vtks_dir) / vtk_file_name;

	pvd_contents += R"(    <DataSet timestep=")" + std::to_string(iteration) + R"(" group="" part="0" file=")"
					+ std::filesystem::relative(file_path, output_dir).string() + R"(" />
)";

	std::ofstream pvd_file(std::filesystem::path(output_dir) / "microenvironment.pvd");

	pvd_file << pvd_contents << R"(  </Collection>
</VTKFile>)";

	pvd_file.close();
}
