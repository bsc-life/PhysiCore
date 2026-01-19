#pragma once

#include <string_view>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <common/generic_agent_solver.h>
#include <common/types.h>
#include <common/vtk_serializer_base.h>

#include "cell.h"
#include "environment.h"

namespace physicore::mechanics::micromechanics {

using common::vtkRealArray;

class vtk_mechanics_serializer : public common::vtk_serializer_base, private generic_agent_solver<agent>
{
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	// Agent Classification
	vtkSmartPointer<vtkRealArray> agent_types_array;
	vtkSmartPointer<vtkRealArray> cell_ids_array;

	// Physics State
	vtkSmartPointer<vtkRealArray> velocities_array;
	vtkSmartPointer<vtkRealArray> previous_velocities_array;
	vtkSmartPointer<vtkRealArray> forces_array;

	// Geometry & Properties
	vtkSmartPointer<vtkRealArray> radii_array;
	vtkSmartPointer<vtkRealArray> is_movable_array;

	// Mechanics Parameters
	vtkSmartPointer<vtkRealArray> cell_cell_adhesion_strength_array;
	vtkSmartPointer<vtkRealArray> cell_cell_repulsion_strength_array;
	vtkSmartPointer<vtkRealArray> relative_maximum_adhesion_distance_array;
	vtkSmartPointer<vtkRealArray> cell_BM_adhesion_strength_array;
	vtkSmartPointer<vtkRealArray> cell_BM_repulsion_strength_array;

	// Mechanics - Attachments
	vtkSmartPointer<vtkRealArray> maximum_number_of_attachments_array;
	vtkSmartPointer<vtkRealArray> attachment_elastic_constant_array;
	vtkSmartPointer<vtkRealArray> attachment_rate_array;
	vtkSmartPointer<vtkRealArray> detachment_rate_array;

	// Mechanics - Advanced (Morse / Kelvin-Voigt)
	vtkSmartPointer<vtkRealArray> cell_residency_array;
	vtkSmartPointer<vtkRealArray> intra_scaling_factors_array;
	vtkSmartPointer<vtkRealArray> intra_equilibrium_distances_array;
	vtkSmartPointer<vtkRealArray> intra_stiffnesses_array;
	vtkSmartPointer<vtkRealArray> spring_constants_array;
	vtkSmartPointer<vtkRealArray> dissipation_rates_array;

	// Motility
	vtkSmartPointer<vtkRealArray> is_motile_array;
	vtkSmartPointer<vtkRealArray> persistence_times_array;
	vtkSmartPointer<vtkRealArray> migration_speeds_array;
	vtkSmartPointer<vtkRealArray> migration_bias_directions_array;
	vtkSmartPointer<vtkRealArray> migration_biases_array;
	vtkSmartPointer<vtkRealArray> motility_directions_array;
	vtkSmartPointer<vtkRealArray> restrict_to_2d_array;
	vtkSmartPointer<vtkRealArray> chemotaxis_index_array;
	vtkSmartPointer<vtkRealArray> chemotaxis_direction_array;

public:
	vtk_mechanics_serializer(std::string_view output_dir);

	void serialize(const environment& e, real_t current_time);
};

} // namespace physicore::mechanics::micromechanics
