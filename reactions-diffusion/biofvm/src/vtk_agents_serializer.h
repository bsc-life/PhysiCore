#pragma once

#include <string_view>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <common/generic_agent_solver.h>
#include <common/types.h>

#include "agent.h"
#include "serializer.h"
#include "vtk_serializer_base.h"

namespace physicore::reactions_diffusion::biofvm {

class vtk_agents_serializer : public vtk_serializer_base, public serializer, private generic_agent_solver<agent>
{
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkSmartPointer<vtkRealArray> volumes_array;
	std::vector<vtkSmartPointer<vtkRealArray>> secretion_rates_arrays;
	std::vector<vtkSmartPointer<vtkRealArray>> saturation_densities_arrays;
	std::vector<vtkSmartPointer<vtkRealArray>> uptake_rates_arrays;
	std::vector<vtkSmartPointer<vtkRealArray>> net_export_rates_arrays;
	std::vector<vtkSmartPointer<vtkRealArray>> internalized_substrates_arrays;
	std::vector<vtkSmartPointer<vtkRealArray>> fraction_released_at_death_arrays;
	std::vector<vtkSmartPointer<vtkRealArray>> fraction_transferred_when_ingested_arrays;

	index_t substrate_count;

public:
	vtk_agents_serializer(std::string_view output_dir, const microenvironment& m);

	void serialize(const microenvironment& m, real_t current_time) override;
};

} // namespace physicore::reactions_diffusion::biofvm
