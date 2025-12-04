#include <filesystem>
#include <iomanip>
#include <sstream>
#include <vtkCellArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

#include <micromechanics/agent_container.h>
#include <micromechanics/environment.h>
#include <micromechanics/vtk_mechanics_serializer.h>

using namespace physicore::mechanics::micromechanics;

vtk_mechanics_serializer::vtk_mechanics_serializer(std::string_view output_dir)
	: common::vtk_serializer_base(output_dir, "vtk_mechanics", "mechanics.pvd")
{
	// Initialize arrays
	velocities_array = vtkSmartPointer<vtkRealArray>::New();
	velocities_array->SetNumberOfComponents(3);
	velocities_array->SetName("velocity");
	unstructured_grid->GetPointData()->AddArray(velocities_array);

	forces_array = vtkSmartPointer<vtkRealArray>::New();
	forces_array->SetNumberOfComponents(3);
	forces_array->SetName("force");
	unstructured_grid->GetPointData()->AddArray(forces_array);

	radii_array = vtkSmartPointer<vtkRealArray>::New();
	radii_array->SetNumberOfComponents(1);
	radii_array->SetName("radius");
	unstructured_grid->GetPointData()->AddArray(radii_array);

	is_movable_array = vtkSmartPointer<vtkRealArray>::New();
	is_movable_array->SetNumberOfComponents(1);
	is_movable_array->SetName("is_movable");
	unstructured_grid->GetPointData()->AddArray(is_movable_array);

	cell_cell_adhesion_strength_array = vtkSmartPointer<vtkRealArray>::New();
	cell_cell_adhesion_strength_array->SetNumberOfComponents(1);
	cell_cell_adhesion_strength_array->SetName("cell_cell_adhesion_strength");
	unstructured_grid->GetPointData()->AddArray(cell_cell_adhesion_strength_array);

	cell_cell_repulsion_strength_array = vtkSmartPointer<vtkRealArray>::New();
	cell_cell_repulsion_strength_array->SetNumberOfComponents(1);
	cell_cell_repulsion_strength_array->SetName("cell_cell_repulsion_strength");
	unstructured_grid->GetPointData()->AddArray(cell_cell_repulsion_strength_array);

	relative_maximum_adhesion_distance_array = vtkSmartPointer<vtkRealArray>::New();
	relative_maximum_adhesion_distance_array->SetNumberOfComponents(1);
	relative_maximum_adhesion_distance_array->SetName("relative_maximum_adhesion_distance");
	unstructured_grid->GetPointData()->AddArray(relative_maximum_adhesion_distance_array);

	is_motile_array = vtkSmartPointer<vtkRealArray>::New();
	is_motile_array->SetNumberOfComponents(1);
	is_motile_array->SetName("is_motile");
	unstructured_grid->GetPointData()->AddArray(is_motile_array);

	motility_vector_array = vtkSmartPointer<vtkRealArray>::New();
	motility_vector_array->SetNumberOfComponents(3);
	motility_vector_array->SetName("motility_vector");
	unstructured_grid->GetPointData()->AddArray(motility_vector_array);

	writer->SetInputData(unstructured_grid);
	writer->SetCompressorTypeToNone();
}

void vtk_mechanics_serializer::serialize(const environment& e, real_t current_time)
{
	const auto& mechanics_data = retrieve_agent_data(*e.agents);
	const auto& base_data = mechanics_data.base_data;

	const index_t agent_count = mechanics_data.agents_count;

	// Create points for agent positions
	auto points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(agent_count);

	// Set point positions from agent data
	for (index_t i = 0; i < agent_count; ++i)
	{
		double pos[3] = { 0.0, 0.0, 0.0 };
		for (index_t d = 0; d < base_data.dims && d < 3; ++d)
		{
			pos[d] = base_data.positions[i * base_data.dims + d];
		}
		points->SetPoint(i, pos);
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
	velocities_array->SetNumberOfTuples(agent_count);
	forces_array->SetNumberOfTuples(agent_count);
	radii_array->SetNumberOfTuples(agent_count);
	is_movable_array->SetNumberOfTuples(agent_count);
	cell_cell_adhesion_strength_array->SetNumberOfTuples(agent_count);
	cell_cell_repulsion_strength_array->SetNumberOfTuples(agent_count);
	relative_maximum_adhesion_distance_array->SetNumberOfTuples(agent_count);
	is_motile_array->SetNumberOfTuples(agent_count);
	motility_vector_array->SetNumberOfTuples(agent_count);

	// Fill in the data
	for (index_t i = 0; i < agent_count; ++i)
	{
		double vel[3] = { 0.0, 0.0, 0.0 };
		double force[3] = { 0.0, 0.0, 0.0 };
		double motility[3] = { 0.0, 0.0, 0.0 };

		for (index_t d = 0; d < base_data.dims && d < 3; ++d)
		{
			vel[d] = mechanics_data.velocities[i * base_data.dims + d];
			force[d] = mechanics_data.forces[i * base_data.dims + d];
			motility[d] = mechanics_data.motility_directions[i * base_data.dims + d];
		}

		velocities_array->SetTuple(i, vel);
		forces_array->SetTuple(i, force);
		motility_vector_array->SetTuple(i, motility);

		radii_array->SetValue(i, mechanics_data.radii[i]);
		is_movable_array->SetValue(i, mechanics_data.is_movable[i]);
		cell_cell_adhesion_strength_array->SetValue(i, mechanics_data.cell_cell_adhesion_strength[i]);
		cell_cell_repulsion_strength_array->SetValue(i, mechanics_data.cell_cell_repulsion_strength[i]);
		relative_maximum_adhesion_distance_array->SetValue(i, mechanics_data.relative_maximum_adhesion_distance[i]);
		is_motile_array->SetValue(i, mechanics_data.is_motile[i]);
	}

	// Write the file
	std::ostringstream ss;
	ss << "mechanics_" << std::setw(6) << std::setfill('0') << iteration << ".vtu";

	auto file_name = ss.str();
	auto file_path = std::filesystem::path(vtks_dir) / file_name;

	writer->SetFileName(file_path.string().c_str());
	writer->Write();

	append_to_pvd(file_name, current_time);

	iteration++;
}
