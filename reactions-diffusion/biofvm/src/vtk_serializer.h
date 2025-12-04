#pragma once

#include <string_view>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

#include <common/types.h>
#include <common/vtk_serializer_base.h>

#include "serializer.h"

namespace physicore::biofvm {

using common::vtkRealArray;

class vtk_serializer : public common::vtk_serializer_base, public serializer
{
	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	vtkSmartPointer<vtkImageData> image_data = vtkSmartPointer<vtkImageData>::New();
	std::vector<vtkSmartPointer<vtkRealArray>> data_arrays;

public:
	vtk_serializer(std::string_view output_dir, microenvironment& m);

	void serialize(const microenvironment& m, real_t current_time) override;
};

} // namespace physicore::biofvm
