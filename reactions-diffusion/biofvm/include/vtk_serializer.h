#pragma once

#include <string_view>
#include <type_traits>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

#include "serializer.h"
#include "types.h"
#include "vtk_serializer_base.h"

namespace physicore::biofvm {

using vtkRealArray = std::conditional_t<std::is_same_v<real_t, float>, vtkFloatArray, vtkDoubleArray>;

class vtk_serializer : public vtk_serializer_base, public serializer
{
	vtkSmartPointer<vtkXMLImageDataWriter> writer;
	vtkSmartPointer<vtkImageData> image_data;
	std::vector<vtkSmartPointer<vtkRealArray>> data_arrays;

public:
	vtk_serializer(std::string_view output_dir, microenvironment& m);

	void serialize(const microenvironment& m, real_t current_time) override;
};

} // namespace physicore::biofvm
