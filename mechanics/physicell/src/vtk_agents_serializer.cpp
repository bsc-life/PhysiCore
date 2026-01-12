#include "vtk_agents_serializer.h"

#include <array>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <utility>
#include <vtkCellArray.h>
#include <vtkCellType.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

using namespace physicore::mechanics::physicell;

vtk_agents_serializer::vtk_agents_serializer(std::string_view output_dir,
											 mechanical_agent_container_interface& container,
											 std::vector<std::string> substrate_names,
											 std::vector<std::string> cell_type_names)
	: vtk_serializer_base(output_dir, "vtk_mechanics_agents", "mechanics_agents.pvd"),
	  container(container),
	  substrate_names(std::move(substrate_names)),
	  cell_type_names(std::move(cell_type_names))
{
	writer->SetInputData(unstructured_grid);
	writer->SetCompressorTypeToNone();
}

std::string vtk_agents_serializer::make_substrate_name(index_t index) const
{
	if (index < substrate_names.size() && !substrate_names[index].empty())
		return substrate_names[index];

	std::ostringstream ss;
	ss << "substrate_" << index;
	return ss.str();
}

std::string vtk_agents_serializer::make_cell_type_name(index_t index) const
{
	if (index < cell_type_names.size() && !cell_type_names[index].empty())
		return cell_type_names[index];

	std::ostringstream ss;
	ss << "cell_type_" << index;
	return ss.str();
}

void vtk_agents_serializer::initialize_arrays(index_t agent_types_count, index_t substrates_count)
{
	radius_array = vtkSmartPointer<vtkRealArray>::New();
	radius_array->SetName("radius");
	radius_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(radius_array);

	cell_cell_adhesion_array = vtkSmartPointer<vtkRealArray>::New();
	cell_cell_adhesion_array->SetName("cell_cell_adhesion_strength");
	cell_cell_adhesion_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(cell_cell_adhesion_array);

	cell_BM_adhesion_array = vtkSmartPointer<vtkRealArray>::New();
	cell_BM_adhesion_array->SetName("cell_BM_adhesion_strength");
	cell_BM_adhesion_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(cell_BM_adhesion_array);

	cell_cell_repulsion_array = vtkSmartPointer<vtkRealArray>::New();
	cell_cell_repulsion_array->SetName("cell_cell_repulsion_strength");
	cell_cell_repulsion_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(cell_cell_repulsion_array);

	cell_BM_repulsion_array = vtkSmartPointer<vtkRealArray>::New();
	cell_BM_repulsion_array->SetName("cell_BM_repulsion_strength");
	cell_BM_repulsion_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(cell_BM_repulsion_array);

	relative_maximum_adhesion_distance_array = vtkSmartPointer<vtkRealArray>::New();
	relative_maximum_adhesion_distance_array->SetName("relative_maximum_adhesion_distance");
	relative_maximum_adhesion_distance_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(relative_maximum_adhesion_distance_array);

	maximum_number_of_attachments_array = vtkSmartPointer<vtkRealArray>::New();
	maximum_number_of_attachments_array->SetName("maximum_number_of_attachments");
	maximum_number_of_attachments_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(maximum_number_of_attachments_array);

	attachment_elastic_constant_array = vtkSmartPointer<vtkRealArray>::New();
	attachment_elastic_constant_array->SetName("attachment_elastic_constant");
	attachment_elastic_constant_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(attachment_elastic_constant_array);

	attachment_rate_array = vtkSmartPointer<vtkRealArray>::New();
	attachment_rate_array->SetName("attachment_rate");
	attachment_rate_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(attachment_rate_array);

	detachment_rate_array = vtkSmartPointer<vtkRealArray>::New();
	detachment_rate_array->SetName("detachment_rate");
	detachment_rate_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(detachment_rate_array);

	is_motile_array = vtkSmartPointer<vtkRealArray>::New();
	is_motile_array->SetName("is_motile");
	is_motile_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(is_motile_array);

	persistence_time_array = vtkSmartPointer<vtkRealArray>::New();
	persistence_time_array->SetName("persistence_time");
	persistence_time_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(persistence_time_array);

	migration_speed_array = vtkSmartPointer<vtkRealArray>::New();
	migration_speed_array->SetName("migration_speed");
	migration_speed_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(migration_speed_array);

	migration_bias_array = vtkSmartPointer<vtkRealArray>::New();
	migration_bias_array->SetName("migration_bias");
	migration_bias_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(migration_bias_array);

	restrict_to_2d_array = vtkSmartPointer<vtkRealArray>::New();
	restrict_to_2d_array->SetName("restrict_to_2d");
	restrict_to_2d_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(restrict_to_2d_array);

	chemotaxis_index_array = vtkSmartPointer<vtkRealArray>::New();
	chemotaxis_index_array->SetName("chemotaxis_index");
	chemotaxis_index_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(chemotaxis_index_array);

	chemotaxis_direction_array = vtkSmartPointer<vtkRealArray>::New();
	chemotaxis_direction_array->SetName("chemotaxis_direction");
	chemotaxis_direction_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(chemotaxis_direction_array);

	simple_pressure_array = vtkSmartPointer<vtkRealArray>::New();
	simple_pressure_array->SetName("simple_pressure");
	simple_pressure_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(simple_pressure_array);

	cell_definition_index_array = vtkSmartPointer<vtkRealArray>::New();
	cell_definition_index_array->SetName("cell_definition_index");
	cell_definition_index_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(cell_definition_index_array);

	is_movable_array = vtkSmartPointer<vtkRealArray>::New();
	is_movable_array->SetName("is_movable");
	is_movable_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(is_movable_array);

	velocity_array = vtkSmartPointer<vtkRealArray>::New();
	velocity_array->SetName("velocity");
	velocity_array->SetNumberOfComponents(3);
	unstructured_grid->GetPointData()->AddArray(velocity_array);

	previous_velocity_array = vtkSmartPointer<vtkRealArray>::New();
	previous_velocity_array->SetName("previous_velocity");
	previous_velocity_array->SetNumberOfComponents(3);
	unstructured_grid->GetPointData()->AddArray(previous_velocity_array);

	migration_bias_direction_array = vtkSmartPointer<vtkRealArray>::New();
	migration_bias_direction_array->SetName("migration_bias_direction");
	migration_bias_direction_array->SetNumberOfComponents(3);
	unstructured_grid->GetPointData()->AddArray(migration_bias_direction_array);

	motility_vector_array = vtkSmartPointer<vtkRealArray>::New();
	motility_vector_array->SetName("motility_vector");
	motility_vector_array->SetNumberOfComponents(3);
	unstructured_grid->GetPointData()->AddArray(motility_vector_array);

	orientation_array = vtkSmartPointer<vtkRealArray>::New();
	orientation_array->SetName("orientation");
	orientation_array->SetNumberOfComponents(3);
	unstructured_grid->GetPointData()->AddArray(orientation_array);

	cell_adhesion_affinity_arrays.clear();
	cell_adhesion_affinity_arrays.reserve(agent_types_count);
	for (index_t i = 0; i < agent_types_count; ++i)
	{
		auto array = vtkSmartPointer<vtkRealArray>::New();
		array->SetNumberOfComponents(1);
		array->SetName((make_cell_type_name(i) + "_cell_adhesion_affinity").c_str());
		cell_adhesion_affinity_arrays.push_back(array);
		unstructured_grid->GetPointData()->AddArray(array);
	}

	chemotactic_sensitivity_arrays.clear();
	chemotactic_sensitivity_arrays.reserve(substrates_count);
	for (index_t i = 0; i < substrates_count; ++i)
	{
		auto array = vtkSmartPointer<vtkRealArray>::New();
		array->SetNumberOfComponents(1);
		array->SetName((make_substrate_name(i) + "_chemotactic_sensitivity").c_str());
		chemotactic_sensitivity_arrays.push_back(array);
		unstructured_grid->GetPointData()->AddArray(array);
	}

	stored_agent_types = agent_types_count;
	stored_substrates = substrates_count;
	initialized = true;
}

void vtk_agents_serializer::serialize(real_t current_time)
{
	auto& data = retrieve_agent_data(container);
	const index_t agent_count = data.agents_count;
	const index_t dims = data.base_data.dims;

	if (!initialized)
	{
		initialize_arrays(data.agent_types_count, data.substrates_count);
	}

	auto points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(static_cast<vtkIdType>(agent_count));

	for (index_t i = 0; i < agent_count; ++i)
	{
		std::array<double, 3> pos = { 0.0, 0.0, 0.0 };
		for (index_t d = 0; d < dims && d < 3; ++d)
		{
			pos[d] = data.base_data.positions[i * dims + d];
		}
		points->SetPoint(static_cast<vtkIdType>(i), pos.data());
	}
	unstructured_grid->SetPoints(points);

	auto cells = vtkSmartPointer<vtkCellArray>::New();
	for (index_t i = 0; i < agent_count; ++i)
	{
		cells->InsertNextCell(1);
		cells->InsertCellPoint(static_cast<vtkIdType>(i));
	}
	unstructured_grid->SetCells(VTK_VERTEX, cells);

	radius_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	cell_cell_adhesion_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	cell_BM_adhesion_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	cell_cell_repulsion_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	cell_BM_repulsion_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	relative_maximum_adhesion_distance_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	maximum_number_of_attachments_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	attachment_elastic_constant_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	attachment_rate_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	detachment_rate_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	is_motile_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	persistence_time_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	migration_speed_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	migration_bias_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	restrict_to_2d_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	chemotaxis_index_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	chemotaxis_direction_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	simple_pressure_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	cell_definition_index_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	is_movable_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));

	velocity_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	previous_velocity_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	migration_bias_direction_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	motility_vector_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	orientation_array->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));

	for (auto& arr : cell_adhesion_affinity_arrays)
	{
		arr->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	}

	for (auto& arr : chemotactic_sensitivity_arrays)
	{
		arr->SetNumberOfTuples(static_cast<vtkIdType>(agent_count));
	}

	for (index_t i = 0; i < agent_count; ++i)
	{
		radius_array->SetValue(static_cast<vtkIdType>(i), data.radius[i]);
		cell_cell_adhesion_array->SetValue(static_cast<vtkIdType>(i),
										   data.mechanics_data.cell_cell_adhesion_strength[i]);
		cell_BM_adhesion_array->SetValue(static_cast<vtkIdType>(i), data.mechanics_data.cell_BM_adhesion_strength[i]);
		cell_cell_repulsion_array->SetValue(static_cast<vtkIdType>(i),
											data.mechanics_data.cell_cell_repulsion_strength[i]);
		cell_BM_repulsion_array->SetValue(static_cast<vtkIdType>(i), data.mechanics_data.cell_BM_repulsion_strength[i]);
		relative_maximum_adhesion_distance_array->SetValue(static_cast<vtkIdType>(i),
														   data.mechanics_data.relative_maximum_adhesion_distance[i]);
		maximum_number_of_attachments_array->SetValue(
			static_cast<vtkIdType>(i), static_cast<real_t>(data.mechanics_data.maximum_number_of_attachments[i]));
		attachment_elastic_constant_array->SetValue(static_cast<vtkIdType>(i),
													data.mechanics_data.attachment_elastic_constant[i]);
		attachment_rate_array->SetValue(static_cast<vtkIdType>(i), data.mechanics_data.attachment_rate[i]);
		detachment_rate_array->SetValue(static_cast<vtkIdType>(i), data.mechanics_data.detachment_rate[i]);

		is_motile_array->SetValue(static_cast<vtkIdType>(i), static_cast<real_t>(data.motility_data.is_motile[i]));
		persistence_time_array->SetValue(static_cast<vtkIdType>(i), data.motility_data.persistence_time[i]);
		migration_speed_array->SetValue(static_cast<vtkIdType>(i), data.motility_data.migration_speed[i]);
		migration_bias_array->SetValue(static_cast<vtkIdType>(i), data.motility_data.migration_bias[i]);
		restrict_to_2d_array->SetValue(static_cast<vtkIdType>(i),
									   static_cast<real_t>(data.motility_data.restrict_to_2d[i]));
		chemotaxis_index_array->SetValue(static_cast<vtkIdType>(i),
										 static_cast<real_t>(data.motility_data.chemotaxis_index[i]));
		chemotaxis_direction_array->SetValue(static_cast<vtkIdType>(i),
											 static_cast<real_t>(data.motility_data.chemotaxis_direction[i]));
		simple_pressure_array->SetValue(static_cast<vtkIdType>(i), data.state_data.simple_pressure[i]);
		cell_definition_index_array->SetValue(static_cast<vtkIdType>(i),
											  static_cast<real_t>(data.state_data.agent_type_index[i]));
		is_movable_array->SetValue(static_cast<vtkIdType>(i), static_cast<real_t>(data.state_data.is_movable[i]));

		std::array<double, 3> velocity = { 0.0, 0.0, 0.0 };
		std::array<double, 3> previous_velocity = { 0.0, 0.0, 0.0 };
		std::array<double, 3> migration_bias_dir = { 0.0, 0.0, 0.0 };
		std::array<double, 3> motility_vector = { 0.0, 0.0, 0.0 };
		std::array<double, 3> orientation = { 0.0, 0.0, 0.0 };

		for (index_t d = 0; d < dims && d < 3; ++d)
		{
			const index_t offset = i * dims + d;
			velocity[d] = data.velocity[offset];
			previous_velocity[d] = data.previous_velocity[offset];
			migration_bias_dir[d] = data.motility_data.migration_bias_direction[offset];
			motility_vector[d] = data.motility_data.motility_vector[offset];
			orientation[d] = data.state_data.orientation[offset];
		}

		velocity_array->SetTuple(static_cast<vtkIdType>(i), velocity.data());
		previous_velocity_array->SetTuple(static_cast<vtkIdType>(i), previous_velocity.data());
		migration_bias_direction_array->SetTuple(static_cast<vtkIdType>(i), migration_bias_dir.data());
		motility_vector_array->SetTuple(static_cast<vtkIdType>(i), motility_vector.data());
		orientation_array->SetTuple(static_cast<vtkIdType>(i), orientation.data());

		for (index_t t = 0; t < stored_agent_types; ++t)
		{
			const index_t idx = i * stored_agent_types + t;
			cell_adhesion_affinity_arrays[t]->SetValue(static_cast<vtkIdType>(i),
													   data.mechanics_data.cell_adhesion_affinities[idx]);
		}

		for (index_t s = 0; s < stored_substrates; ++s)
		{
			const index_t idx = i * stored_substrates + s;
			chemotactic_sensitivity_arrays[s]->SetValue(static_cast<vtkIdType>(i),
														data.motility_data.chemotactic_sensitivities[idx]);
		}
	}

	std::ostringstream ss;
	ss << "mechanics_agents_" << std::setw(6) << std::setfill('0') << iteration << ".vtu";

	auto file_name = ss.str();
	auto file_path = std::filesystem::path(vtks_dir) / file_name;

	writer->SetFileName(file_path.string().c_str());
	writer->Write();

	append_to_pvd(file_name, current_time);

	iteration++;
}
