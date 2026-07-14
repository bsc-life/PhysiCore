#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <pugixml.hpp>

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

		const std::array<sindex_t, 3> bounding_box_mins = { 0, 0, 0 };
		const std::array<sindex_t, 3> bounding_box_maxs = { (sindex_t)grid_shape[0] * (sindex_t)voxel_shape[0],
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

	// Parse all values from a named CellData DataArray in a .vti file.
	static std::vector<real_t> read_cell_array(const std::filesystem::path& vti_path, const char* array_name)
	{
		pugi::xml_document doc;
		if (!doc.load_file(vti_path.string().c_str())) return {};
		auto cell_data = doc.child("VTKFile").child("ImageData").child("Piece").child("CellData");
		for (auto da : cell_data.children("DataArray"))
			if (std::string(da.attribute("Name").value()) == array_name)
			{
				std::vector<real_t> vals;
				std::istringstream ss(da.text().get());
				real_t v;
				while (ss >> v) vals.push_back(v);
				return vals;
			}
		return {};
	}
};

TEST_F(VtkSerializerTest, ConstructorInitialization)
{
	auto m = create_test_microenvironment();

	// Test constructor doesn't throw
	EXPECT_NO_THROW({ const vtk_serializer serializer(test_output_dir.string(), *m); });

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
	const std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

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

	// Read and validate VTK file structure using pugixml
	auto vtk_dir = test_output_dir / "vtk_microenvironment";
	auto vti_file = vtk_dir / "microenvironment_000000.vti";

	pugi::xml_document doc;
	ASSERT_TRUE(doc.load_file(vti_file.string().c_str()));

	auto image_data_node = doc.child("VTKFile").child("ImageData");

	// Check dimensions from WholeExtent (dims = extent[2i+1] - extent[2i] + 1)
	std::string extent_str = image_data_node.attribute("WholeExtent").value();
	int e0, e1, e2, e3, e4, e5;
	std::istringstream(extent_str) >> e0 >> e1 >> e2 >> e3 >> e4 >> e5;
	EXPECT_EQ(e1 - e0 + 1, 4);
	EXPECT_EQ(e3 - e2 + 1, 4);
	EXPECT_EQ(e5 - e4 + 1, 4);

	// Check spacing
	std::string spacing_str = image_data_node.attribute("Spacing").value();
	double sx, sy, sz;
	std::istringstream(spacing_str) >> sx >> sy >> sz;
	EXPECT_DOUBLE_EQ(sx, 20.0);
	EXPECT_DOUBLE_EQ(sy, 20.0);
	EXPECT_DOUBLE_EQ(sz, 20.0);

	// Check CellData arrays exist and have correct metadata
	auto cell_data = image_data_node.child("Piece").child("CellData");
	ASSERT_FALSE(cell_data.empty());

	auto find_array = [&](const char* name) -> pugi::xml_node {
		for (auto da : cell_data.children("DataArray"))
			if (std::string(da.attribute("Name").value()) == name) return da;
		return {};
	};

	auto o2_node = find_array("O2");
	ASSERT_FALSE(o2_node.empty()) << "O2 array not found";
	EXPECT_EQ(std::stoi(o2_node.attribute("NumberOfComponents").value()), 1);
	EXPECT_EQ(std::stoull(o2_node.attribute("NumberOfTuples").value()), 27u);

	auto glucose_node = find_array("Glucose");
	ASSERT_FALSE(glucose_node.empty()) << "Glucose array not found";
	EXPECT_EQ(std::stoi(glucose_node.attribute("NumberOfComponents").value()), 1);
	EXPECT_EQ(std::stoull(glucose_node.attribute("NumberOfTuples").value()), 27u);
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

	// Verify single substrate array via pugixml
	auto vtk_dir = test_output_dir / "vtk_microenvironment";
	auto vti_file = vtk_dir / "microenvironment_000000.vti";

	pugi::xml_document doc;
	ASSERT_TRUE(doc.load_file(vti_file.string().c_str()));
	auto cell_data = doc.child("VTKFile").child("ImageData").child("Piece").child("CellData");

	int n_arrays = 0;
	for (auto da : cell_data.children("DataArray")) { ++n_arrays; (void)da; }
	EXPECT_EQ(n_arrays, 1);

	bool found_o2 = false;
	for (auto da : cell_data.children("DataArray"))
		if (std::string(da.attribute("Name").value()) == "O2") { found_o2 = true; break; }
	EXPECT_TRUE(found_o2);
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
		builder.add_density(name, "unit", 1.0, 0.01, static_cast<real_t>(10 + i)); // Different initial conditions

		// Add boundary conditions for each substrate
		builder.add_boundary_dirichlet_conditions(i, // substrate index
												  { static_cast<real_t>(15 + i), static_cast<real_t>(14 + i),
													static_cast<real_t>(13 + i) }, // min boundary values
												  { static_cast<real_t>(20 + i), static_cast<real_t>(19 + i),
													static_cast<real_t>(18 + i) }, // max boundary values
												  { true, true, true },			   // min boundary conditions
												  { true, true, true }			   // max boundary conditions
		);
	}

	auto m = builder.build();
	m->solver->initialize(*m); // Initialize densities

	vtk_serializer serializer(test_output_dir.string(), *m);

	EXPECT_NO_THROW(serializer.serialize(*m, 0.0));

	// Verify all substrate arrays via pugixml
	auto vtk_dir = test_output_dir / "vtk_microenvironment";
	auto vti_file = vtk_dir / "microenvironment_000000.vti";

	pugi::xml_document doc;
	ASSERT_TRUE(doc.load_file(vti_file.string().c_str()));
	auto cell_data = doc.child("VTKFile").child("ImageData").child("Piece").child("CellData");

	int n_arrays = 0;
	for (auto da : cell_data.children("DataArray")) { ++n_arrays; (void)da; }
	EXPECT_EQ(n_arrays, (int)substrate_names.size());

	for (const auto& name : substrate_names)
	{
		bool found = false;
		for (auto da : cell_data.children("DataArray"))
			if (std::string(da.attribute("Name").value()) == name) { found = true; break; }
		EXPECT_TRUE(found) << "Substrate array '" << name << "' not found";
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

	// Verify extent calculation with non-zero bounding box mins via pugixml
	auto vtk_dir = test_output_dir / "vtk_microenvironment";
	auto vti_file = vtk_dir / "microenvironment_000000.vti";

	pugi::xml_document doc;
	ASSERT_TRUE(doc.load_file(vti_file.string().c_str()));
	auto image_data_node = doc.child("VTKFile").child("ImageData");

	std::string extent_str = image_data_node.attribute("WholeExtent").value();
	int e0, e1, e2, e3, e4, e5;
	std::istringstream(extent_str) >> e0 >> e1 >> e2 >> e3 >> e4 >> e5;

	// Check extent calculation: bounding_box_mins / voxel_shape
	EXPECT_EQ(e0, 5);  // 100/20
	EXPECT_EQ(e1, 8);  // 5 + 3 grid cells
	EXPECT_EQ(e2, 10); // 200/20
	EXPECT_EQ(e3, 13); // 10 + 3 grid cells
	EXPECT_EQ(e4, 15); // 300/20
	EXPECT_EQ(e5, 18); // 15 + 3 grid cells
}

TEST_F(VtkSerializerTest, VtkRealArrayTypeConsistency)
{
	auto m = create_test_microenvironment();
	vtk_serializer serializer(test_output_dir.string(), *m);

	serializer.serialize(*m, 0.0);

	auto vtk_dir = test_output_dir / "vtk_microenvironment";
	auto vti_file = vtk_dir / "microenvironment_000000.vti";

	pugi::xml_document doc;
	ASSERT_TRUE(doc.load_file(vti_file.string().c_str()));
	auto cell_data = doc.child("VTKFile").child("ImageData").child("Piece").child("CellData");

	for (auto da : cell_data.children("DataArray"))
	{
		if (std::string(da.attribute("Name").value()) == "O2")
		{
			const char* type_str = da.attribute("type").value();
			if (std::is_same_v<real_t, float>)
				EXPECT_STREQ(type_str, "Float32");
			else
				EXPECT_STREQ(type_str, "Float64");
			break;
		}
	}
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

	// Read back and check specific voxel values via pugixml
	auto o2_vals = read_cell_array(vti_file, "O2");
	ASSERT_EQ(o2_vals.size(), m->mesh.voxel_count());

	auto glucose_vals = read_cell_array(vti_file, "Glucose");
	ASSERT_EQ(glucose_vals.size(), m->mesh.voxel_count());

	for (index_t z = 0; z < m->mesh.grid_shape[2]; ++z)
		for (index_t y = 0; y < m->mesh.grid_shape[1]; ++y)
			for (index_t x = 0; x < m->mesh.grid_shape[0]; ++x)
			{
				const std::size_t voxel_idx = m->mesh.linearize(x, y, z);
				const real_t value = o2_vals[voxel_idx];

				if (x == 0)
					EXPECT_EQ(value, 100.0); // High boundary
				else if (x == m->mesh.grid_shape[0] - 1)
					EXPECT_EQ(value, 5.0); // Low boundary
				else
					EXPECT_EQ(value, 20.0);

				EXPECT_EQ(glucose_vals[voxel_idx], m->get_substrate_density(1, x, y, z));
			}

	m->solver->solve(*m, 1);
	serializer.serialize(*m, 0.1);

	auto vti_file2 = vtk_dir / "microenvironment_000001.vti";
	EXPECT_TRUE(std::filesystem::exists(vti_file2));

	auto glucose_vals2 = read_cell_array(vti_file2, "Glucose");
	ASSERT_EQ(glucose_vals2.size(), m->mesh.voxel_count());

	for (index_t z = 0; z < m->mesh.grid_shape[2]; ++z)
		for (index_t y = 0; y < m->mesh.grid_shape[1]; ++y)
			for (index_t x = 0; x < m->mesh.grid_shape[0]; ++x)
			{
				const std::size_t voxel_idx = m->mesh.linearize(x, y, z);
				EXPECT_EQ(glucose_vals2[voxel_idx], m->get_substrate_density(1, x, y, z));
			}
}
