#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>
#include <vtkXMLImageDataReader.h>

#include <gtest/gtest.h>

#include "microenvironment.h"
#include "microenvironment_builder.h"
#include "vtk_serializer.h"

using namespace physicore;
using namespace physicore::biofvm;

class VtkSerializerTest : public ::testing::Test
{
protected:
	void SetUp() override
	{
		// Create a temporary output directory for tests
		test_output_dir = std::filesystem::temp_directory_path() / "vtk_serializer_test";
		std::filesystem::create_directories(test_output_dir);

		// Clean up any existing files
		if (std::filesystem::exists(test_output_dir))
		{
			std::filesystem::remove_all(test_output_dir);
		}
		std::filesystem::create_directories(test_output_dir);
	}

	void TearDown() override
	{
		// Clean up test files
		if (std::filesystem::exists(test_output_dir))
		{
			std::filesystem::remove_all(test_output_dir);
		}
	}

	static std::unique_ptr<microenvironment> create_test_microenvironment(
		index_t dims = 3, std::array<index_t, 3> grid_shape = { 3, 3, 3 },
		std::array<index_t, 3> voxel_shape = { 20, 20, 20 })
	{
		microenvironment_builder builder;
		builder.set_name("test_env");
		builder.set_time_units("min");
		builder.set_space_units("um");
		builder.set_time_step(0.01);

		std::array<sindex_t, 3> bounding_box_mins = { 0, 0, 0 };
		std::array<sindex_t, 3> bounding_box_maxs = { (sindex_t)grid_shape[0] * (sindex_t)voxel_shape[0],
													  (sindex_t)grid_shape[1] * (sindex_t)voxel_shape[1],
													  (sindex_t)grid_shape[2] * (sindex_t)voxel_shape[2] };

		builder.resize(dims, bounding_box_mins, bounding_box_maxs, voxel_shape);

		// Add substrates with specific initial conditions
		builder.add_density("O2", "mmHg", 1.0, 0.01, 38.0);	  // Typical oxygen concentration
		builder.add_density("Glucose", "mM", 0.5, 0.02, 5.5); // Typical glucose concentration

		auto m = builder.build();

		// Run solver to initialize substrate densities
		m->solver->solve(*m, 1);

		return m;
	}

	std::filesystem::path test_output_dir;
};

TEST_F(VtkSerializerTest, ConstructorInitialization)
{
	auto m = create_test_microenvironment();

	// Test constructor doesn't throw
	EXPECT_NO_THROW({ vtk_serializer serializer(test_output_dir.string(), *m); });

	// Check that directories are created
	auto vtk_dir = test_output_dir / "vtk_microenvironment";
	EXPECT_TRUE(std::filesystem::exists(vtk_dir));
}

TEST_F(VtkSerializerTest, SerializeCreatesFiles)
{
	auto m = create_test_microenvironment();
	vtk_serializer serializer(test_output_dir.string(), *m);

	// Serialize once
	EXPECT_NO_THROW(serializer.serialize(*m, 0.0));

	// Check that VTK file is created
	auto vtk_dir = test_output_dir / "vtk_microenvironment";
	auto vti_file = vtk_dir / "microenvironment_000000.vti";
	EXPECT_TRUE(std::filesystem::exists(vti_file));

	// Check that PVD file is created
	auto pvd_file = test_output_dir / "microenvironment.pvd";
	EXPECT_TRUE(std::filesystem::exists(pvd_file));
}

TEST_F(VtkSerializerTest, SerializeMultipleTimes)
{
	auto m = create_test_microenvironment();
	vtk_serializer serializer(test_output_dir.string(), *m);

	// Serialize multiple times
	for (int i = 0; i < 3; ++i)
	{
		EXPECT_NO_THROW(serializer.serialize(*m, 0.0));
	}

	auto vtk_dir = test_output_dir / "vtk_microenvironment";

	// Check that multiple VTK files are created
	EXPECT_TRUE(std::filesystem::exists(vtk_dir / "microenvironment_000000.vti"));
	EXPECT_TRUE(std::filesystem::exists(vtk_dir / "microenvironment_000001.vti"));
	EXPECT_TRUE(std::filesystem::exists(vtk_dir / "microenvironment_000002.vti"));
}

TEST_F(VtkSerializerTest, PvdFileContainsCorrectEntries)
{
	auto m = create_test_microenvironment();
	vtk_serializer serializer(test_output_dir.string(), *m);

	// Serialize twice
	serializer.serialize(*m, 0.1);
	serializer.serialize(*m, 0.2);

	// Read PVD file content
	auto pvd_file = test_output_dir / "microenvironment.pvd";
	ASSERT_TRUE(std::filesystem::exists(pvd_file));

	std::ifstream file(pvd_file);
	std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

	// Check XML structure
	EXPECT_TRUE(content.find("<?xml version=\"1.0\"?>") != std::string::npos);
	EXPECT_TRUE(content.find("<VTKFile type=\"Collection\"") != std::string::npos);
	EXPECT_TRUE(content.find("<Collection>") != std::string::npos);
	EXPECT_TRUE(content.find("</Collection>") != std::string::npos);
	EXPECT_TRUE(content.find("</VTKFile>") != std::string::npos);

	// Check timestep entries
	EXPECT_TRUE(content.find("timestep=\"0.1") != std::string::npos);
	EXPECT_TRUE(content.find("timestep=\"0.2") != std::string::npos);

	// Check file references
	EXPECT_TRUE(content.find("microenvironment_000000.vti") != std::string::npos);
	EXPECT_TRUE(content.find("microenvironment_000001.vti") != std::string::npos);
}

TEST_F(VtkSerializerTest, VtkFileStructure)
{
	auto m = create_test_microenvironment();
	vtk_serializer serializer(test_output_dir.string(), *m);

	serializer.serialize(*m, 0.0);

	// Read and validate VTK file structure using VTK reader
	auto vtk_dir = test_output_dir / "vtk_microenvironment";
	auto vti_file = vtk_dir / "microenvironment_000000.vti";

	auto reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	reader->SetFileName(vti_file.string().c_str());
	reader->Update();

	auto* image_data = reader->GetOutput();
	ASSERT_NE(image_data, nullptr);

	// Check dimensions
	std::array<int, 3> dims {};
	image_data->GetDimensions(dims.data());
	EXPECT_EQ(dims[0], 4); // grid_shape[0] + 1
	EXPECT_EQ(dims[1], 4); // grid_shape[1] + 1
	EXPECT_EQ(dims[2], 4); // grid_shape[2] + 1

	// Check spacing
	std::array<double, 3> spacing {};
	image_data->GetSpacing(spacing.data());
	EXPECT_DOUBLE_EQ(spacing[0], 20.0);
	EXPECT_DOUBLE_EQ(spacing[1], 20.0);
	EXPECT_DOUBLE_EQ(spacing[2], 20.0);

	// Check that arrays exist for substrates
	auto* cell_data = image_data->GetCellData();
	ASSERT_NE(cell_data, nullptr);

	EXPECT_NE(cell_data->GetArray("O2"), nullptr);
	EXPECT_NE(cell_data->GetArray("Glucose"), nullptr);

	// Check array properties
	auto* o2_array = cell_data->GetArray("O2");
	EXPECT_EQ(o2_array->GetNumberOfComponents(), 1);
	EXPECT_EQ(o2_array->GetNumberOfTuples(), 27); // 3x3x3 = 27 voxels

	auto* glucose_array = cell_data->GetArray("Glucose");
	EXPECT_EQ(glucose_array->GetNumberOfComponents(), 1);
	EXPECT_EQ(glucose_array->GetNumberOfTuples(), 27); // 3x3x3 = 27 voxels
}

TEST_F(VtkSerializerTest, HandleDifferentMeshDimensions)
{
	// Test 1D mesh
	{
		auto m1d = create_test_microenvironment(1, { 5, 1, 1 }, { 10, 1, 1 });
		EXPECT_NO_THROW({
			vtk_serializer serializer(test_output_dir.string(), *m1d);
			serializer.serialize(*m1d, 0.0);
		});
	}

	// Test 2D mesh
	{
		auto m2d = create_test_microenvironment(2, { 4, 4, 1 }, { 15, 15, 1 });
		EXPECT_NO_THROW({
			vtk_serializer serializer(test_output_dir.string(), *m2d);
			serializer.serialize(*m2d, 0.0);
		});
	}
}

TEST_F(VtkSerializerTest, SingleSubstrate)
{
	microenvironment_builder builder;
	builder.set_name("single_substrate_env");
	builder.resize(3, { 0, 0, 0 }, { 40, 40, 40 }, { 20, 20, 20 });
	builder.add_density("O2", "mmHg", 1.0, 0.01, 38.0);

	// Add boundary conditions for single substrate
	builder.add_boundary_dirichlet_conditions(0,					// O2 density index
											  { 40.0, 40.0, 40.0 }, // min boundary values
											  { 40.0, 40.0, 40.0 }, // max boundary values
											  { true, true, true }, // min boundary conditions
											  { true, true, true }	// max boundary conditions
	);

	auto m = builder.build();
	m->solver->initialize(*m); // Initialize densities

	vtk_serializer serializer(test_output_dir.string(), *m);

	EXPECT_NO_THROW(serializer.serialize(*m, 0.0));

	// Verify single substrate array
	auto vtk_dir = test_output_dir / "vtk_microenvironment";
	auto vti_file = vtk_dir / "microenvironment_000000.vti";

	auto reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	reader->SetFileName(vti_file.string().c_str());
	reader->Update();

	auto* image_data = reader->GetOutput();
	auto* cell_data = image_data->GetCellData();

	EXPECT_EQ(cell_data->GetNumberOfArrays(), 1);
	EXPECT_NE(cell_data->GetArray("O2"), nullptr);
}

TEST_F(VtkSerializerTest, ManySubstrates)
{
	microenvironment_builder builder;
	builder.set_name("many_substrates_env");
	builder.resize(3, { 0, 0, 0 }, { 40, 40, 40 }, { 20, 20, 20 });

	// Add many substrates
	std::vector<std::string> substrate_names = { "O2", "Glucose", "Lactate", "ATP", "CO2", "H2O" };

	for (size_t i = 0; i < substrate_names.size(); ++i)
	{
		const auto& name = substrate_names[i];
		builder.add_density(name, "unit", 1.0, 0.01, 10.0 + i); // Different initial conditions

		// Add boundary conditions for each substrate
		builder.add_boundary_dirichlet_conditions(i, // substrate index
												  { static_cast<real_t>(15.0 + i), static_cast<real_t>(14.0 + i),
													static_cast<real_t>(13.0 + i) }, // min boundary values
												  { static_cast<real_t>(20.0 + i), static_cast<real_t>(19.0 + i),
													static_cast<real_t>(18.0 + i) }, // max boundary values
												  { true, true, true },				 // min boundary conditions
												  { true, true, true }				 // max boundary conditions
		);
	}

	auto m = builder.build();
	m->solver->initialize(*m); // Initialize densities

	vtk_serializer serializer(test_output_dir.string(), *m);

	EXPECT_NO_THROW(serializer.serialize(*m, 0.0));

	// Verify all substrate arrays
	auto vtk_dir = test_output_dir / "vtk_microenvironment";
	auto vti_file = vtk_dir / "microenvironment_000000.vti";

	auto reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	reader->SetFileName(vti_file.string().c_str());
	reader->Update();

	auto* image_data = reader->GetOutput();
	auto* cell_data = image_data->GetCellData();

	EXPECT_EQ(cell_data->GetNumberOfArrays(), substrate_names.size());

	for (const auto& name : substrate_names)
	{
		EXPECT_NE(cell_data->GetArray(name.c_str()), nullptr) << "Substrate array '" << name << "' not found";
	}
}

TEST_F(VtkSerializerTest, NonZeroBoundingBoxMins)
{
	microenvironment_builder builder;
	builder.set_name("offset_env");
	builder.resize(3, { 100, 200, 300 }, { 160, 260, 360 }, { 20, 20, 20 });
	builder.add_density("O2", "mmHg", 1.0, 0.01, 35.0);

	// Add boundary conditions with offset coordinates
	builder.add_boundary_dirichlet_conditions(0,					// O2 density index
											  { 38.0, 36.0, 34.0 }, // min boundary values
											  { 42.0, 40.0, 38.0 }, // max boundary values
											  { true, true, true }, // min boundary conditions
											  { true, true, true }	// max boundary conditions
	);

	auto m = builder.build();
	m->solver->initialize(*m); // Initialize densities

	vtk_serializer serializer(test_output_dir.string(), *m);

	EXPECT_NO_THROW(serializer.serialize(*m, 0.0));

	// Verify extent calculation with non-zero bounding box mins
	auto vtk_dir = test_output_dir / "vtk_microenvironment";
	auto vti_file = vtk_dir / "microenvironment_000000.vti";

	auto reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	reader->SetFileName(vti_file.string().c_str());
	reader->Update();

	auto* image_data = reader->GetOutput();

	std::array<int, 6> extent {};
	image_data->GetExtent(extent.data());

	// Check extent calculation: bounding_box_mins / voxel_shape
	EXPECT_EQ(extent[0], 5);  // 100/20
	EXPECT_EQ(extent[1], 8);  // 5 + 3 grid cells
	EXPECT_EQ(extent[2], 10); // 200/20
	EXPECT_EQ(extent[3], 13); // 10 + 3 grid cells
	EXPECT_EQ(extent[4], 15); // 300/20
	EXPECT_EQ(extent[5], 18); // 15 + 3 grid cells
}

TEST_F(VtkSerializerTest, VtkRealArrayTypeConsistency)
{
	auto m = create_test_microenvironment();
	vtk_serializer serializer(test_output_dir.string(), *m);

	serializer.serialize(*m, 0.0);

	auto vtk_dir = test_output_dir / "vtk_microenvironment";
	auto vti_file = vtk_dir / "microenvironment_000000.vti";

	auto reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	reader->SetFileName(vti_file.string().c_str());
	reader->Update();

	auto* image_data = reader->GetOutput();
	auto* cell_data = image_data->GetCellData();
	auto* array = cell_data->GetArray("O2");

	// Check that array type matches real_t (double in this case)
	if (std::is_same_v<real_t, float>)
		EXPECT_EQ(array->GetDataType(), VTK_FLOAT);
	else
		EXPECT_EQ(array->GetDataType(), VTK_DOUBLE);
}

TEST_F(VtkSerializerTest, BoundaryConditionsEffect)
{
	// Create a microenvironment with strong boundary conditions to test serialization of gradients
	microenvironment_builder builder;
	builder.set_name("boundary_test_env");
	builder.resize(3, { 0, 0, 0 }, { 60, 60, 60 }, { 20, 20, 20 });
	builder.add_density("O2", "mmHg", 0, 0, 20.0);		  // Zero diffusion and decay
	builder.add_density("Glucose", "mM", 0.5, 0.02, 5.5); // Typical glucose concentration

	// Set high concentration on one boundary, low on the opposite
	builder.add_boundary_dirichlet_conditions(0,					  // O2 density index
											  { 100.0, 20.0, 20.0 },  // high x_min, normal y_min, z_min
											  { 5.0, 20.0, 20.0 },	  // low x_max, normal y_max, z_max
											  { true, false, false }, // only x_min boundary active
											  { true, false, false }  // only x_max boundary active
	);

	auto m = builder.build();

	m->solver->solve(*m, 1);

	vtk_serializer serializer(test_output_dir.string(), *m);
	EXPECT_NO_THROW(serializer.serialize(*m, 0.0));

	// Verify files exist
	auto vtk_dir = test_output_dir / "vtk_microenvironment";
	auto vti_file = vtk_dir / "microenvironment_000000.vti";
	EXPECT_TRUE(std::filesystem::exists(vti_file));

	// Read back the VTI file and check specific voxel values
	auto reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	reader->SetFileName(vti_file.string().c_str());
	reader->Update();

	auto* image_data = reader->GetOutput();
	ASSERT_NE(image_data, nullptr);

	auto* cell_data = image_data->GetCellData();
	ASSERT_NE(cell_data, nullptr);

	auto* o2_array = cell_data->GetArray("O2");
	ASSERT_NE(o2_array, nullptr);

	EXPECT_EQ(o2_array->GetNumberOfTuples(), m->mesh.voxel_count());

	for (index_t z = 0; z < m->mesh.grid_shape[2]; ++z)
		for (index_t y = 0; y < m->mesh.grid_shape[1]; ++y)
			for (index_t x = 0; x < m->mesh.grid_shape[0]; ++x)
			{
				std::size_t voxel_idx = m->mesh.linearize(x, y, z);
				real_t value = o2_array->GetTuple1(voxel_idx);

				if (x == 0)
					EXPECT_EQ(value, 100.0); // High boundary
				else if (x == m->mesh.grid_shape[0] - 1)
					EXPECT_EQ(value, 5.0); // Low boundary
				else
					EXPECT_EQ(value, 20.0);

				real_t value2 = cell_data->GetArray("Glucose")->GetTuple1(voxel_idx);
				EXPECT_EQ(value2, m->get_substrate_density(1, x, y, z));
			}

	m->solver->solve(*m, 1);
	serializer.serialize(*m, 0.1);

	auto vti_file2 = vtk_dir / "microenvironment_000001.vti";
	EXPECT_TRUE(std::filesystem::exists(vti_file2));

	// Read back the VTI file and check specific voxel values
	reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	reader->SetFileName(vti_file2.string().c_str());
	reader->Update();

	image_data = reader->GetOutput();
	ASSERT_NE(image_data, nullptr);

	cell_data = image_data->GetCellData();
	ASSERT_NE(cell_data, nullptr);

	auto* glucose_array = cell_data->GetArray("Glucose");
	ASSERT_NE(glucose_array, nullptr);

	EXPECT_EQ(glucose_array->GetNumberOfTuples(), m->mesh.voxel_count());

	for (index_t z = 0; z < m->mesh.grid_shape[2]; ++z)
		for (index_t y = 0; y < m->mesh.grid_shape[1]; ++y)
			for (index_t x = 0; x < m->mesh.grid_shape[0]; ++x)
			{
				std::size_t voxel_idx = m->mesh.linearize(x, y, z);
				real_t value = glucose_array->GetTuple1(voxel_idx);
				EXPECT_EQ(value, m->get_substrate_density(1, x, y, z));
			}
}
