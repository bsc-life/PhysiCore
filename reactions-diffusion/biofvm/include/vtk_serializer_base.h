#pragma once

#include <string>
#include <string_view>

namespace physicore::biofvm {

class vtk_serializer_base
{
protected:
	std::size_t iteration;
	std::string output_dir;
	std::string vtks_dir;

	std::string pvd_contents;

	void append_to_pvd(std::string_view vtk_file_name);

public:
	vtk_serializer_base(std::string_view output_dir, std::string_view vtks_dir_name);
};

} // namespace physicore::biofvm
