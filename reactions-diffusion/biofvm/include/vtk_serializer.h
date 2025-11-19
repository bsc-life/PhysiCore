#pragma once

#include <string_view>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

#include "serializer.h"
#include "types.h"
#include "vtk_serializer_base.h"

namespace physicore::biofvm {

class vtk_serializer : public vtk_serializer_base, public serializer
{
	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	vtkSmartPointer<vtkImageData> image_data = vtkSmartPointer<vtkImageData>::New();
	std::vector<vtkSmartPointer<vtkRealArray>> data_arrays;

public:
	vtk_serializer(std::string_view output_dir, microenvironment& m);

	void serialize(const microenvironment& m, real_t current_time) override;
};

} // namespace physicore::biofvm
