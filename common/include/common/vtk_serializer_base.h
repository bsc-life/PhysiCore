#pragma once

#include <string>
#include <string_view>
#include <type_traits>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>

#include <common/types.h>

namespace physicore::common {

using vtkRealArray = std::conditional_t<std::is_same_v<real_t, float>, vtkFloatArray, vtkDoubleArray>;

class vtk_serializer_base
{
	std::size_t iteration_ = 0;
	std::string output_dir_;
	std::string vtks_dir_;
	std::string pvd_file_name_;

	std::string pvd_contents_;

protected:
	[[nodiscard]] std::size_t iteration() const { return iteration_; }
	void advance_iteration() { ++iteration_; }

	[[nodiscard]] const std::string& vtks_dir() const { return vtks_dir_; }

	void append_to_pvd(std::string_view vtk_file_name, real_t current_time);

public:
	vtk_serializer_base(std::string_view output_dir, std::string_view vtks_dir_name, std::string_view pvd_file_name);
};

} // namespace physicore::common
