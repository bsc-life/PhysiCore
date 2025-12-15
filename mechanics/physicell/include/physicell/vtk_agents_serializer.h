#pragma once

#include <string>
#include <string_view>
#include <vector>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <common/generic_agent_solver.h>
#include <common/types.h>

#include "mechanical_agent.h"
#include "mechanical_agent_container.h"
#include "serializer.h"
#include "vtk_serializer_base.h"

namespace physicore::mechanics::physicell {

class vtk_agents_serializer : public vtk_serializer_base,
							  public serializer,
							  private generic_agent_solver<mechanical_agent>
{
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	mechanical_agent_container_interface& container;

	vtkSmartPointer<vtkRealArray> radius_array;
	vtkSmartPointer<vtkRealArray> cell_cell_adhesion_array;
	vtkSmartPointer<vtkRealArray> cell_BM_adhesion_array;
	vtkSmartPointer<vtkRealArray> cell_cell_repulsion_array;
	vtkSmartPointer<vtkRealArray> cell_BM_repulsion_array;
	vtkSmartPointer<vtkRealArray> relative_maximum_adhesion_distance_array;
	vtkSmartPointer<vtkRealArray> maximum_number_of_attachments_array;
	vtkSmartPointer<vtkRealArray> attachment_elastic_constant_array;
	vtkSmartPointer<vtkRealArray> attachment_rate_array;
	vtkSmartPointer<vtkRealArray> detachment_rate_array;

	vtkSmartPointer<vtkRealArray> is_motile_array;
	vtkSmartPointer<vtkRealArray> persistence_time_array;
	vtkSmartPointer<vtkRealArray> migration_speed_array;
	vtkSmartPointer<vtkRealArray> migration_bias_array;
	vtkSmartPointer<vtkRealArray> restrict_to_2d_array;
	vtkSmartPointer<vtkRealArray> chemotaxis_index_array;
	vtkSmartPointer<vtkRealArray> chemotaxis_direction_array;
	vtkSmartPointer<vtkRealArray> simple_pressure_array;
	vtkSmartPointer<vtkRealArray> cell_definition_index_array;
	vtkSmartPointer<vtkRealArray> is_movable_array;

	vtkSmartPointer<vtkRealArray> velocity_array;
	vtkSmartPointer<vtkRealArray> previous_velocity_array;
	vtkSmartPointer<vtkRealArray> migration_bias_direction_array;
	vtkSmartPointer<vtkRealArray> motility_vector_array;
	vtkSmartPointer<vtkRealArray> orientation_array;

	std::vector<vtkSmartPointer<vtkRealArray>> cell_adhesion_affinity_arrays;
	std::vector<vtkSmartPointer<vtkRealArray>> chemotactic_sensitivity_arrays;

	std::vector<std::string> substrate_names;
	std::vector<std::string> cell_type_names;

	index_t stored_agent_types = 0;
	index_t stored_substrates = 0;
	bool initialized = false;

	void initialize_arrays(index_t agent_types_count, index_t substrates_count);
	std::string make_substrate_name(index_t index) const;
	std::string make_cell_type_name(index_t index) const;

public:
	vtk_agents_serializer(std::string_view output_dir, mechanical_agent_container_interface& container,
						  std::vector<std::string> substrate_names = {},
						  std::vector<std::string> cell_type_names = {});

	void serialize(real_t current_time) override;
};

} // namespace physicore::mechanics::physicell
