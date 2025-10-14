#include "vtk_serializer.h"

#include <filesystem>
#include <sstream>
#include <vtkCellData.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#include "microenvironment.h"

using namespace physicore::biofvm;

vtk_serializer::vtk_serializer(std::string_view output_dir, microenvironment& m)
	: vtk_serializer_base(output_dir, "vtk_microenvironment")
{
	auto x_extent_start = m.mesh.bounding_box_mins[0] / (sindex_t)m.mesh.voxel_shape[0];
	auto x_extent_end = x_extent_start + (sindex_t)m.mesh.grid_shape[0];

	auto y_extent_start = m.mesh.bounding_box_mins[1] / (sindex_t)m.mesh.voxel_shape[1];
	auto y_extent_end = y_extent_start + (sindex_t)m.mesh.grid_shape[1];

	auto z_extent_start = m.mesh.bounding_box_mins[2] / (sindex_t)m.mesh.voxel_shape[2];
	auto z_extent_end = z_extent_start + (sindex_t)m.mesh.grid_shape[2];

	if (m.mesh.dims == 3)
	{
		image_data->SetExtent(x_extent_start, x_extent_end, y_extent_start, y_extent_end, z_extent_start, z_extent_end);
		image_data->SetSpacing(m.mesh.voxel_shape[0], m.mesh.voxel_shape[1], m.mesh.voxel_shape[2]);
	}
	else
	{
		image_data->SetExtent(x_extent_start, x_extent_end, y_extent_start, y_extent_end, 0, 0);
		image_data->SetSpacing(m.mesh.voxel_shape[0], m.mesh.voxel_shape[1], 0);
	}

	data_arrays.reserve(m.substrates_count);
	assert(m.substrates_names.size() == m.substrates_count);
	for (index_t i = 0; i < m.substrates_count; ++i)
	{
		data_arrays.emplace_back(vtkSmartPointer<vtkRealArray>::New());
		data_arrays[i]->SetNumberOfComponents(1);
		data_arrays[i]->SetNumberOfTuples(m.mesh.voxel_count());
		data_arrays[i]->SetName(m.substrates_names[i].c_str());
		image_data->GetCellData()->AddArray(data_arrays[i]);
	}

	writer->SetInputData(image_data);
	writer->SetCompressorTypeToNone();
}

void vtk_serializer::serialize(const microenvironment& m, real_t current_time)
{
	for (index_t z_idx = 0; z_idx < m.mesh.grid_shape[2]; ++z_idx)
		for (index_t y_idx = 0; y_idx < m.mesh.grid_shape[1]; ++y_idx)
			for (index_t x_idx = 0; x_idx < m.mesh.grid_shape[0]; ++x_idx)
				for (index_t s_idx = 0; s_idx < m.substrates_count; ++s_idx)
				{
					std::size_t voxel_idx = m.mesh.linearize(x_idx, y_idx, z_idx);
					data_arrays[s_idx]->SetValue(voxel_idx, m.get_substrate_density(s_idx, x_idx, y_idx, z_idx));
				}

	std::ostringstream ss;

	ss << "microenvironment_" << std::setw(6) << std::setfill('0') << iteration << ".vti";

	auto file_name = ss.str();
	auto file_path = std::filesystem::path(vtks_dir) / file_name;

	writer->SetFileName(file_path.string().c_str());
	writer->Write();

	append_to_pvd(file_name, current_time);

	iteration++;
}
