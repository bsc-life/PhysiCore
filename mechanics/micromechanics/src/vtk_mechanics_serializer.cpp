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
	agent_types_array = vtkSmartPointer<vtkRealArray>::New();
	agent_types_array->SetNumberOfComponents(1);
	agent_types_array->SetName("agent_type");
	unstructured_grid->GetPointData()->AddArray(agent_types_array);

	cell_ids_array = vtkSmartPointer<vtkRealArray>::New();
	cell_ids_array->SetNumberOfComponents(1);
	cell_ids_array->SetName("cell_id");
	unstructured_grid->GetPointData()->AddArray(cell_ids_array);

	// Initialize arrays
	velocities_array = vtkSmartPointer<vtkRealArray>::New();
	velocities_array->SetNumberOfComponents(3);
	velocities_array->SetName("velocity");
	unstructured_grid->GetPointData()->AddArray(velocities_array);

	previous_velocities_array = vtkSmartPointer<vtkRealArray>::New();
	previous_velocities_array->SetNumberOfComponents(3);
	previous_velocities_array->SetName("previous_velocity");
	unstructured_grid->GetPointData()->AddArray(previous_velocities_array);

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

	cell_BM_adhesion_strength_array = vtkSmartPointer<vtkRealArray>::New();
	cell_BM_adhesion_strength_array->SetNumberOfComponents(1);
	cell_BM_adhesion_strength_array->SetName("cell_BM_adhesion_strength");
	unstructured_grid->GetPointData()->AddArray(cell_BM_adhesion_strength_array);

	cell_BM_repulsion_strength_array = vtkSmartPointer<vtkRealArray>::New();
	cell_BM_repulsion_strength_array->SetNumberOfComponents(1);
	cell_BM_repulsion_strength_array->SetName("cell_BM_repulsion_strength");
	unstructured_grid->GetPointData()->AddArray(cell_BM_repulsion_strength_array);

	maximum_number_of_attachments_array = vtkSmartPointer<vtkRealArray>::New();
	maximum_number_of_attachments_array->SetNumberOfComponents(1);
	maximum_number_of_attachments_array->SetName("maximum_number_of_attachments");
	unstructured_grid->GetPointData()->AddArray(maximum_number_of_attachments_array);

	attachment_elastic_constant_array = vtkSmartPointer<vtkRealArray>::New();
	attachment_elastic_constant_array->SetNumberOfComponents(1);
	attachment_elastic_constant_array->SetName("attachment_elastic_constant");
	unstructured_grid->GetPointData()->AddArray(attachment_elastic_constant_array);

	attachment_rate_array = vtkSmartPointer<vtkRealArray>::New();
	attachment_rate_array->SetNumberOfComponents(1);
	attachment_rate_array->SetName("attachment_rate");
	unstructured_grid->GetPointData()->AddArray(attachment_rate_array);

	detachment_rate_array = vtkSmartPointer<vtkRealArray>::New();
	detachment_rate_array->SetNumberOfComponents(1);
	detachment_rate_array->SetName("detachment_rate");
	unstructured_grid->GetPointData()->AddArray(detachment_rate_array);

	cell_residency_array = vtkSmartPointer<vtkRealArray>::New();
	cell_residency_array->SetNumberOfComponents(1);
	cell_residency_array->SetName("cell_residency");
	unstructured_grid->GetPointData()->AddArray(cell_residency_array);

	intra_scaling_factors_array = vtkSmartPointer<vtkRealArray>::New();
	intra_scaling_factors_array->SetNumberOfComponents(1);
	intra_scaling_factors_array->SetName("intra_scaling_factor");
	unstructured_grid->GetPointData()->AddArray(intra_scaling_factors_array);

	intra_equilibrium_distances_array = vtkSmartPointer<vtkRealArray>::New();
	intra_equilibrium_distances_array->SetNumberOfComponents(1);
	intra_equilibrium_distances_array->SetName("intra_equilibrium_distance");
	unstructured_grid->GetPointData()->AddArray(intra_equilibrium_distances_array);

	intra_stiffnesses_array = vtkSmartPointer<vtkRealArray>::New();
	intra_stiffnesses_array->SetNumberOfComponents(1);
	intra_stiffnesses_array->SetName("intra_stiffness");
	unstructured_grid->GetPointData()->AddArray(intra_stiffnesses_array);

	spring_constants_array = vtkSmartPointer<vtkRealArray>::New();
	spring_constants_array->SetNumberOfComponents(1);
	spring_constants_array->SetName("spring_constant");
	unstructured_grid->GetPointData()->AddArray(spring_constants_array);

	dissipation_rates_array = vtkSmartPointer<vtkRealArray>::New();
	dissipation_rates_array->SetNumberOfComponents(1);
	dissipation_rates_array->SetName("dissipation_rate");
	unstructured_grid->GetPointData()->AddArray(dissipation_rates_array);

	is_motile_array = vtkSmartPointer<vtkRealArray>::New();
	is_motile_array->SetNumberOfComponents(1);
	is_motile_array->SetName("is_motile");
	unstructured_grid->GetPointData()->AddArray(is_motile_array);

	persistence_times_array = vtkSmartPointer<vtkRealArray>::New();
	persistence_times_array->SetNumberOfComponents(1);
	persistence_times_array->SetName("persistence_time");
	unstructured_grid->GetPointData()->AddArray(persistence_times_array);

	migration_speeds_array = vtkSmartPointer<vtkRealArray>::New();
	migration_speeds_array->SetNumberOfComponents(1);
	migration_speeds_array->SetName("migration_speed");
	unstructured_grid->GetPointData()->AddArray(migration_speeds_array);

	migration_bias_directions_array = vtkSmartPointer<vtkRealArray>::New();
	migration_bias_directions_array->SetNumberOfComponents(3);
	migration_bias_directions_array->SetName("migration_bias_direction");
	unstructured_grid->GetPointData()->AddArray(migration_bias_directions_array);

	migration_biases_array = vtkSmartPointer<vtkRealArray>::New();
	migration_biases_array->SetNumberOfComponents(1);
	migration_biases_array->SetName("migration_bias");
	unstructured_grid->GetPointData()->AddArray(migration_biases_array);

	motility_directions_array = vtkSmartPointer<vtkRealArray>::New();
	motility_directions_array->SetNumberOfComponents(3);
	motility_directions_array->SetName("motility_direction");
	unstructured_grid->GetPointData()->AddArray(motility_directions_array);

	restrict_to_2d_array = vtkSmartPointer<vtkRealArray>::New();
	restrict_to_2d_array->SetNumberOfComponents(1);
	restrict_to_2d_array->SetName("restrict_to_2d");
	unstructured_grid->GetPointData()->AddArray(restrict_to_2d_array);

	chemotaxis_index_array = vtkSmartPointer<vtkRealArray>::New();
	chemotaxis_index_array->SetNumberOfComponents(1);
	chemotaxis_index_array->SetName("chemotaxis_index");
	unstructured_grid->GetPointData()->AddArray(chemotaxis_index_array);

	chemotaxis_direction_array = vtkSmartPointer<vtkRealArray>::New();
	chemotaxis_direction_array->SetNumberOfComponents(1);
	chemotaxis_direction_array->SetName("chemotaxis_direction");
	unstructured_grid->GetPointData()->AddArray(chemotaxis_direction_array);

	writer->SetInputData(unstructured_grid);
	writer->SetCompressorTypeToNone();
}

void vtk_mechanics_serializer::serialize(const environment& e, real_t current_time)
{
	const auto& mechanics_data = retrieve_agent_data(*e.agents);
	const auto& base_data = mechanics_data.base_data;

	const index_t agent_count = mechanics_data.agents_count;
	const auto vtk_agent_count = static_cast<vtkIdType>(agent_count);

	// Create points for agent positions
	auto points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(vtk_agent_count);

	// Set point positions from agent data
	for (index_t i = 0; i < agent_count; ++i)
	{
		std::array<double, 3> pos = { 0.0, 0.0, 0.0 };
		for (index_t d = 0; d < base_data.dims && d < 3; ++d)
		{
			pos[d] = base_data.positions[i * base_data.dims + d];
		}
		points->SetPoint(static_cast<vtkIdType>(i), pos.data());
	}
	unstructured_grid->SetPoints(points);

	// Create vertex cells (one per agent)
	auto cells = vtkSmartPointer<vtkCellArray>::New();
	for (index_t i = 0; i < agent_count; ++i)
	{
		cells->InsertNextCell(1);
		cells->InsertCellPoint(static_cast<vtkIdType>(i));
	}
	unstructured_grid->SetCells(VTK_VERTEX, cells);

	// Set array sizes
	agent_types_array->SetNumberOfTuples(vtk_agent_count);
	cell_ids_array->SetNumberOfTuples(vtk_agent_count);
	velocities_array->SetNumberOfTuples(vtk_agent_count);
	previous_velocities_array->SetNumberOfTuples(vtk_agent_count);
	forces_array->SetNumberOfTuples(vtk_agent_count);
	radii_array->SetNumberOfTuples(vtk_agent_count);
	is_movable_array->SetNumberOfTuples(vtk_agent_count);
	cell_cell_adhesion_strength_array->SetNumberOfTuples(vtk_agent_count);
	cell_cell_repulsion_strength_array->SetNumberOfTuples(vtk_agent_count);
	relative_maximum_adhesion_distance_array->SetNumberOfTuples(vtk_agent_count);
	cell_BM_adhesion_strength_array->SetNumberOfTuples(vtk_agent_count);
	cell_BM_repulsion_strength_array->SetNumberOfTuples(vtk_agent_count);
	maximum_number_of_attachments_array->SetNumberOfTuples(vtk_agent_count);
	attachment_elastic_constant_array->SetNumberOfTuples(vtk_agent_count);
	attachment_rate_array->SetNumberOfTuples(vtk_agent_count);
	detachment_rate_array->SetNumberOfTuples(vtk_agent_count);
	cell_residency_array->SetNumberOfTuples(vtk_agent_count);
	intra_scaling_factors_array->SetNumberOfTuples(vtk_agent_count);
	intra_equilibrium_distances_array->SetNumberOfTuples(vtk_agent_count);
	intra_stiffnesses_array->SetNumberOfTuples(vtk_agent_count);
	spring_constants_array->SetNumberOfTuples(vtk_agent_count);
	dissipation_rates_array->SetNumberOfTuples(vtk_agent_count);
	is_motile_array->SetNumberOfTuples(vtk_agent_count);
	persistence_times_array->SetNumberOfTuples(vtk_agent_count);
	migration_speeds_array->SetNumberOfTuples(vtk_agent_count);
	migration_bias_directions_array->SetNumberOfTuples(vtk_agent_count);
	migration_biases_array->SetNumberOfTuples(vtk_agent_count);
	motility_directions_array->SetNumberOfTuples(vtk_agent_count);
	restrict_to_2d_array->SetNumberOfTuples(vtk_agent_count);
	chemotaxis_index_array->SetNumberOfTuples(vtk_agent_count);
	chemotaxis_direction_array->SetNumberOfTuples(vtk_agent_count);

	// Fill in the data
	for (index_t i = 0; i < agent_count; ++i)
	{
		std::array<double, 3> vel = { 0.0, 0.0, 0.0 };
		std::array<double, 3> prev_vel = { 0.0, 0.0, 0.0 };
		std::array<double, 3> force = { 0.0, 0.0, 0.0 };
		std::array<double, 3> migration_bias_dir = { 0.0, 0.0, 0.0 };
		std::array<double, 3> motility_dir = { 0.0, 0.0, 0.0 };

		for (index_t d = 0; d < base_data.dims && d < 3; ++d)
		{
			vel[d] = mechanics_data.velocities[i * base_data.dims + d];
			prev_vel[d] = mechanics_data.previous_velocities[i * base_data.dims + d];
			force[d] = mechanics_data.forces[i * base_data.dims + d];
			migration_bias_dir[d] = mechanics_data.migration_bias_directions[i * base_data.dims + d];
			motility_dir[d] = mechanics_data.motility_directions[i * base_data.dims + d];
		}

		agent_types_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.agent_types[i]);
		cell_ids_array->SetValue(static_cast<vtkIdType>(i), static_cast<double>(mechanics_data.cell_ids[i]));

		velocities_array->SetTuple(static_cast<vtkIdType>(i), vel.data());
		previous_velocities_array->SetTuple(static_cast<vtkIdType>(i), prev_vel.data());
		forces_array->SetTuple(static_cast<vtkIdType>(i), force.data());
		migration_bias_directions_array->SetTuple(static_cast<vtkIdType>(i), migration_bias_dir.data());
		motility_directions_array->SetTuple(static_cast<vtkIdType>(i), motility_dir.data());

		radii_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.radii[i]);
		is_movable_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.is_movable[i]);
		cell_cell_adhesion_strength_array->SetValue(static_cast<vtkIdType>(i),
											   mechanics_data.cell_cell_adhesion_strength[i]);
		cell_cell_repulsion_strength_array->SetValue(static_cast<vtkIdType>(i),
												mechanics_data.cell_cell_repulsion_strength[i]);
		relative_maximum_adhesion_distance_array->SetValue(static_cast<vtkIdType>(i),
													 mechanics_data.relative_maximum_adhesion_distance[i]);
		cell_BM_adhesion_strength_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.cell_BM_adhesion_strength[i]);
		cell_BM_repulsion_strength_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.cell_BM_repulsion_strength[i]);
		maximum_number_of_attachments_array->SetValue(static_cast<vtkIdType>(i),
													mechanics_data.maximum_number_of_attachments[i]);
		attachment_elastic_constant_array->SetValue(static_cast<vtkIdType>(i),
													 mechanics_data.attachment_elastic_constant[i]);
		attachment_rate_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.attachment_rate[i]);
		detachment_rate_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.detachment_rate[i]);
		cell_residency_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.cell_residency[i]);
		intra_scaling_factors_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.intra_scaling_factors[i]);
		intra_equilibrium_distances_array->SetValue(static_cast<vtkIdType>(i),
													mechanics_data.intra_equilibrium_distances[i]);
		intra_stiffnesses_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.intra_stiffnesses[i]);
		spring_constants_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.spring_constants[i]);
		dissipation_rates_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.dissipation_rates[i]);
		is_motile_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.is_motile[i]);
		persistence_times_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.persistence_times[i]);
		migration_speeds_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.migration_speeds[i]);
		migration_biases_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.migration_biases[i]);
		restrict_to_2d_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.restrict_to_2d[i]);
		chemotaxis_index_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.chemotaxis_index[i]);
		chemotaxis_direction_array->SetValue(static_cast<vtkIdType>(i), mechanics_data.chemotaxis_direction[i]);
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
