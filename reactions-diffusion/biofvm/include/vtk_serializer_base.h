#pragma once

#include <string>
#include <string_view>
#include <type_traits>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>

#include "types.h"

namespace physicore::biofvm {

using vtkRealArray = std::conditional_t<std::is_same_v<real_t, float>, vtkFloatArray, vtkDoubleArray>;

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
