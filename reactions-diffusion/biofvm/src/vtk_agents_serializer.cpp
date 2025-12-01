#include "vtk_agents_serializer.h"

#include <filesystem>
#include <iomanip>
#include <sstream>
#include <vtkCellArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

#include "agent_container.h"
#include "microenvironment.h"

using namespace physicore::biofvm;

vtk_agents_serializer::vtk_agents_serializer(std::string_view output_dir, const microenvironment& m)
	: vtk_serializer_base(output_dir, "vtk_agents", "agents.pvd"), substrate_count(m.substrates_count)
{
	// Initialize volumes array
	volumes_array = vtkSmartPointer<vtkRealArray>::New();
	volumes_array->SetNumberOfComponents(1);
	volumes_array->SetName("volume");
	unstructured_grid->GetPointData()->AddArray(volumes_array);

	// Initialize substrate-related arrays (one array per substrate)
	for (index_t i = 0; i < substrate_count; ++i)
	{
		auto secretion_array = vtkSmartPointer<vtkRealArray>::New();
		secretion_array->SetNumberOfComponents(1);
		secretion_array->SetName((m.substrates_names[i] + "_secretion_rate").c_str());
		secretion_rates_arrays.push_back(secretion_array);
		unstructured_grid->GetPointData()->AddArray(secretion_array);

		auto saturation_array = vtkSmartPointer<vtkRealArray>::New();
		saturation_array->SetNumberOfComponents(1);
		saturation_array->SetName((m.substrates_names[i] + "_saturation_density").c_str());
		saturation_densities_arrays.push_back(saturation_array);
		unstructured_grid->GetPointData()->AddArray(saturation_array);

		auto uptake_array = vtkSmartPointer<vtkRealArray>::New();
		uptake_array->SetNumberOfComponents(1);
		uptake_array->SetName((m.substrates_names[i] + "_uptake_rate").c_str());
		uptake_rates_arrays.push_back(uptake_array);
		unstructured_grid->GetPointData()->AddArray(uptake_array);

		auto net_export_array = vtkSmartPointer<vtkRealArray>::New();
		net_export_array->SetNumberOfComponents(1);
		net_export_array->SetName((m.substrates_names[i] + "_net_export_rate").c_str());
		net_export_rates_arrays.push_back(net_export_array);
		unstructured_grid->GetPointData()->AddArray(net_export_array);

		auto internalized_array = vtkSmartPointer<vtkRealArray>::New();
		internalized_array->SetNumberOfComponents(1);
		internalized_array->SetName((m.substrates_names[i] + "_internalized_substrate").c_str());
		internalized_substrates_arrays.push_back(internalized_array);
		unstructured_grid->GetPointData()->AddArray(internalized_array);

		auto fraction_released_array = vtkSmartPointer<vtkRealArray>::New();
		fraction_released_array->SetNumberOfComponents(1);
		fraction_released_array->SetName((m.substrates_names[i] + "_fraction_released_at_death").c_str());
		fraction_released_at_death_arrays.push_back(fraction_released_array);
		unstructured_grid->GetPointData()->AddArray(fraction_released_array);

		auto fraction_transferred_array = vtkSmartPointer<vtkRealArray>::New();
		fraction_transferred_array->SetNumberOfComponents(1);
		fraction_transferred_array->SetName((m.substrates_names[i] + "_fraction_transferred_when_ingested").c_str());
		fraction_transferred_when_ingested_arrays.push_back(fraction_transferred_array);
		unstructured_grid->GetPointData()->AddArray(fraction_transferred_array);
	}

	writer->SetInputData(unstructured_grid);
	writer->SetCompressorTypeToNone();
}

void vtk_agents_serializer::serialize(const microenvironment& m, real_t current_time)
{
	const auto& biofvm_data = retrieve_agent_data(*m.agents);
	const auto& base_data = biofvm_data.base_data;

	const index_t agent_count = biofvm_data.agents_count;

	// Create points for agent positions
	auto points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(agent_count);

	// Set point positions from agent data
	for (index_t i = 0; i < agent_count; ++i)
	{
		std::array<double, 3> pos = { 0.0, 0.0, 0.0 };
		for (index_t d = 0; d < base_data.dims && d < 3; ++d)
		{
			pos[d] = base_data.positions[i * base_data.dims + d];
		}
		points->SetPoint(i, pos.data());
	}
	unstructured_grid->SetPoints(points);

	// Create vertex cells (one per agent)
	auto cells = vtkSmartPointer<vtkCellArray>::New();
	for (index_t i = 0; i < agent_count; ++i)
	{
		cells->InsertNextCell(1);
		cells->InsertCellPoint(i);
	}
	unstructured_grid->SetCells(VTK_VERTEX, cells);

	// Set array sizes
	volumes_array->SetNumberOfTuples(agent_count);
	for (index_t s = 0; s < substrate_count; ++s)
	{
		secretion_rates_arrays[s]->SetNumberOfTuples(agent_count);
		saturation_densities_arrays[s]->SetNumberOfTuples(agent_count);
		uptake_rates_arrays[s]->SetNumberOfTuples(agent_count);
		net_export_rates_arrays[s]->SetNumberOfTuples(agent_count);
		internalized_substrates_arrays[s]->SetNumberOfTuples(agent_count);
		fraction_released_at_death_arrays[s]->SetNumberOfTuples(agent_count);
		fraction_transferred_when_ingested_arrays[s]->SetNumberOfTuples(agent_count);
	}

	// Fill in the data
	for (index_t i = 0; i < agent_count; ++i)
	{
		volumes_array->SetValue(i, biofvm_data.volumes[i]);

		for (index_t s = 0; s < substrate_count; ++s)
		{
			index_t idx = i * substrate_count + s;
			secretion_rates_arrays[s]->SetValue(i, biofvm_data.secretion_rates[idx]);
			saturation_densities_arrays[s]->SetValue(i, biofvm_data.saturation_densities[idx]);
			uptake_rates_arrays[s]->SetValue(i, biofvm_data.uptake_rates[idx]);
			net_export_rates_arrays[s]->SetValue(i, biofvm_data.net_export_rates[idx]);
			internalized_substrates_arrays[s]->SetValue(i, biofvm_data.internalized_substrates[idx]);
			fraction_released_at_death_arrays[s]->SetValue(i, biofvm_data.fraction_released_at_death[idx]);
			fraction_transferred_when_ingested_arrays[s]->SetValue(i,
																   biofvm_data.fraction_transferred_when_ingested[idx]);
		}
	}

	// Write the file
	std::ostringstream ss;
	ss << "agents_" << std::setw(6) << std::setfill('0') << iteration << ".vtu";

	auto file_name = ss.str();
	auto file_path = std::filesystem::path(vtks_dir) / file_name;

	writer->SetFileName(file_path.string().c_str());
	writer->Write();

	append_to_pvd(file_name, current_time);

	iteration++;
}
