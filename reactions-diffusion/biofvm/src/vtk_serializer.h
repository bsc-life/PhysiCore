#pragma once

#include <array>
#include <string>
#include <string_view>
#include <vector>

#include <common/types.h>

#include "serializer.h"
#include "vtk_serializer_base.h"

namespace physicore::biofvm {

class vtk_serializer : public vtk_serializer_base, public serializer
{
	std::array<int, 6> extent_; // xmin xmax ymin ymax zmin zmax
	std::array<double, 3> spacing_;
	index_t substrates_count_;
	std::vector<std::string> substrate_names_;

public:
	vtk_serializer(std::string_view output_dir, microenvironment& m);

	void serialize(const microenvironment& m, real_t current_time) override;
};

} // namespace physicore::biofvm
