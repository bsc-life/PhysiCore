#pragma once

#include <string>
#include <string_view>

#include "types.h"

namespace physicore::biofvm {

class vtk_serializer_base
{
protected:
	std::size_t iteration = 0;
	std::string output_dir;
	std::string vtks_dir;
	std::string pvd_file_name;

	std::string pvd_contents;

	void append_to_pvd(std::string_view vtk_file_name, real_t current_time);

public:
	vtk_serializer_base(std::string_view output_dir, std::string_view vtks_dir_name, std::string_view pvd_file_name);
};

} // namespace physicore::biofvm
