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

	// Physics State
	vtkSmartPointer<vtkRealArray> velocities_array;
	vtkSmartPointer<vtkRealArray> forces_array;

	// Geometry & Properties
	vtkSmartPointer<vtkRealArray> radii_array;
	vtkSmartPointer<vtkRealArray> is_movable_array;

	// Mechanics Parameters
	vtkSmartPointer<vtkRealArray> cell_cell_adhesion_strength_array;
	vtkSmartPointer<vtkRealArray> cell_cell_repulsion_strength_array;
	vtkSmartPointer<vtkRealArray> relative_maximum_adhesion_distance_array;

	// Motility
	vtkSmartPointer<vtkRealArray> is_motile_array;
	vtkSmartPointer<vtkRealArray> motility_vector_array;

public:
	vtk_mechanics_serializer(std::string_view output_dir);

	void serialize(const environment& e, real_t current_time);
};

} // namespace physicore::mechanics::micromechanics
