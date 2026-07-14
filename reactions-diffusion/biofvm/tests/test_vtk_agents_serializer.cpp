#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <pugixml.hpp>

#include <gtest/gtest.h>

#include "agent_container.h"
#include "microenvironment.h"
#include "microenvironment_builder.h"
#include "vtk_agents_serializer.h"

using namespace physicore;
using namespace physicore::biofvm;

class VtkAgentsSerializerTest : public ::testing::Test
{
protected:
	void SetUp() override
	{
		// Create a temporary output directory for tests
		test_output_dir = std::filesystem::temp_directory_path() / "vtk_agents_serializer_test";
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

	static std::unique_ptr<microenvironment> create_test_microenvironment()
	{
		microenvironment_builder builder;
		builder.set_name("test_env");
		builder.set_time_units("min");
		builder.set_space_units("um");
		builder.set_time_step(0.01);

		const std::array<sindex_t, 3> bounding_box_mins = { 0, 0, 0 };
		const std::array<sindex_t, 3> bounding_box_maxs = { 60, 60, 60 };
		const std::array<index_t, 3> voxel_shape = { 20, 20, 20 };

		builder.resize(3, bounding_box_mins, bounding_box_maxs, voxel_shape);

		// Add substrates
		builder.add_density("O2", "mmHg", 1.0, 0.01, 38.0);
		builder.add_density("Glucose", "mM", 0.5, 0.02, 5.5);

		return builder.build();
	}

	std::filesystem::path test_output_dir;

	// Parse all values from a named PointData DataArray in a .vtu file.
	static std::vector<real_t> read_point_array(const std::filesystem::path& vtu_path, const char* array_name)
	{
		pugi::xml_document doc;
		if (!doc.load_file(vtu_path.string().c_str())) return {};
		auto point_data =
			doc.child("VTKFile").child("UnstructuredGrid").child("Piece").child("PointData");
		for (auto da : point_data.children("DataArray"))
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

	// Parse agent positions (3-component Points DataArray) from a .vtu file.
	static std::vector<std::array<double, 3>> read_points(const std::filesystem::path& vtu_path)
	{
		pugi::xml_document doc;
		if (!doc.load_file(vtu_path.string().c_str())) return {};
		auto points_da =
			doc.child("VTKFile").child("UnstructuredGrid").child("Piece").child("Points").child("DataArray");
		std::vector<std::array<double, 3>> pts;
		std::istringstream ss(points_da.text().get());
		double x, y, z;
		while (ss >> x >> y >> z) pts.push_back({ x, y, z });
		return pts;
	}
};

TEST_F(VtkAgentsSerializerTest, ConstructorInitialization)
{
	auto m = create_test_microenvironment();

	// Test constructor doesn't throw
	EXPECT_NO_THROW({ const vtk_agents_serializer serializer(test_output_dir.string(), *m); });

	// Check that directories are created
	auto vtk_dir = test_output_dir / "vtk_agents";
	EXPECT_TRUE(std::filesystem::exists(vtk_dir));
}

TEST_F(VtkAgentsSerializerTest, SerializeWithNoAgents)
{
	auto m = create_test_microenvironment();
	vtk_agents_serializer serializer(test_output_dir.string(), *m);

	// Serialize with no agents
	EXPECT_NO_THROW(serializer.serialize(*m, 0.0));

	// Check that VTK file is created
	auto vtk_dir = test_output_dir / "vtk_agents";
	auto vtu_file = vtk_dir / "agents_000000.vtu";
	EXPECT_TRUE(std::filesystem::exists(vtu_file));

	// Check that PVD file is created
	auto pvd_file = test_output_dir / "agents.pvd";
	EXPECT_TRUE(std::filesystem::exists(pvd_file));
}

TEST_F(VtkAgentsSerializerTest, SerializeWithSingleAgent)
{
	auto m = create_test_microenvironment();

	// Create an agent
	m->agents->create();
	auto container = std::dynamic_pointer_cast<agent_container>(m->agents);
	ASSERT_NE(container, nullptr);

	auto& base_data = std::get<0>(container->agent_datas);
	auto& biofvm_data = std::get<1>(container->agent_datas);

	// Set agent position
	base_data->positions[0] = 10.0;
	base_data->positions[1] = 20.0;
	base_data->positions[2] = 30.0;

	// Set agent volume
	biofvm_data->volumes[0] = 100.0;

	// Set agent substrate data
	biofvm_data->secretion_rates[0] = 1.0; // O2
	biofvm_data->secretion_rates[1] = 2.0; // Glucose
	biofvm_data->saturation_densities[0] = 3.0;
	biofvm_data->saturation_densities[1] = 4.0;

	vtk_agents_serializer serializer(test_output_dir.string(), *m);
	serializer.serialize(*m, 0.0);

	// Read and validate VTK file structure via pugixml
	auto vtk_dir = test_output_dir / "vtk_agents";
	auto vtu_file = vtk_dir / "agents_000000.vtu";

	pugi::xml_document doc;
	ASSERT_TRUE(doc.load_file(vtu_file.string().c_str()));
	auto piece = doc.child("VTKFile").child("UnstructuredGrid").child("Piece");

	// Check number of points
	EXPECT_EQ(std::stoi(piece.attribute("NumberOfPoints").value()), 1);

	// Check point position
	auto pts = read_points(vtu_file);
	ASSERT_EQ(pts.size(), 1u);
	EXPECT_DOUBLE_EQ(pts[0][0], 10.0);
	EXPECT_DOUBLE_EQ(pts[0][1], 20.0);
	EXPECT_DOUBLE_EQ(pts[0][2], 30.0);

	// Check PointData arrays
	auto vol = read_point_array(vtu_file, "volume");
	ASSERT_EQ(vol.size(), 1u);
	EXPECT_DOUBLE_EQ(static_cast<double>(vol[0]), 100.0);

	auto o2_sec = read_point_array(vtu_file, "O2_secretion_rate");
	ASSERT_EQ(o2_sec.size(), 1u);
	EXPECT_DOUBLE_EQ(static_cast<double>(o2_sec[0]), 1.0);

	auto glc_sec = read_point_array(vtu_file, "Glucose_secretion_rate");
	ASSERT_EQ(glc_sec.size(), 1u);
	EXPECT_DOUBLE_EQ(static_cast<double>(glc_sec[0]), 2.0);

	auto o2_sat = read_point_array(vtu_file, "O2_saturation_density");
	ASSERT_EQ(o2_sat.size(), 1u);
	EXPECT_DOUBLE_EQ(static_cast<double>(o2_sat[0]), 3.0);

	auto glc_sat = read_point_array(vtu_file, "Glucose_saturation_density");
	ASSERT_EQ(glc_sat.size(), 1u);
	EXPECT_DOUBLE_EQ(static_cast<double>(glc_sat[0]), 4.0);
}

TEST_F(VtkAgentsSerializerTest, SerializeWithMultipleAgents)
{
	auto m = create_test_microenvironment();

	// Create multiple agents
	const int num_agents = 5;
	for (int i = 0; i < num_agents; ++i)
	{
		m->agents->create();
	}

	auto container = std::dynamic_pointer_cast<agent_container>(m->agents);
	ASSERT_NE(container, nullptr);

	auto& base_data = std::get<0>(container->agent_datas);
	auto& biofvm_data = std::get<1>(container->agent_datas);

	// Set agent data
	for (int i = 0; i < num_agents; ++i)
	{
		base_data->positions[i * 3 + 0] = i * 10.0;
		base_data->positions[i * 3 + 1] = i * 20.0;
		base_data->positions[i * 3 + 2] = i * 30.0;

		biofvm_data->volumes[i] = (i + 1) * 100.0;

		// Set substrate data for each agent
		for (index_t s = 0; s < m->substrates_count; ++s)
		{
			const index_t idx = i * m->substrates_count + s;
			biofvm_data->secretion_rates[idx] = static_cast<real_t>(i + s + 1);
			biofvm_data->saturation_densities[idx] = static_cast<real_t>((i + 1) * (s + 1) * 10);
			biofvm_data->uptake_rates[idx] = (i + 1) * 0.5;
			biofvm_data->net_export_rates[idx] = (i + 1) * 0.1;
			biofvm_data->internalized_substrates[idx] = i * 5.0;
			biofvm_data->fraction_released_at_death[idx] = 0.5 + i * 0.1;
			biofvm_data->fraction_transferred_when_ingested[idx] = 0.3 + i * 0.05;
		}
	}

	vtk_agents_serializer serializer(test_output_dir.string(), *m);
	serializer.serialize(*m, 0.0);

	// Read and validate VTK file via pugixml
	auto vtk_dir = test_output_dir / "vtk_agents";
	auto vtu_file = vtk_dir / "agents_000000.vtu";

	pugi::xml_document doc;
	ASSERT_TRUE(doc.load_file(vtu_file.string().c_str()));
	auto piece = doc.child("VTKFile").child("UnstructuredGrid").child("Piece");

	// Check number of points
	EXPECT_EQ(std::stoi(piece.attribute("NumberOfPoints").value()), num_agents);

	// Check point positions
	auto pts = read_points(vtu_file);
	ASSERT_EQ((int)pts.size(), num_agents);
	for (int i = 0; i < num_agents; ++i)
	{
		EXPECT_DOUBLE_EQ(pts[i][0], i * 10.0);
		EXPECT_DOUBLE_EQ(pts[i][1], i * 20.0);
		EXPECT_DOUBLE_EQ(pts[i][2], i * 30.0);
	}

	// Check per-agent data arrays
	auto vol = read_point_array(vtu_file, "volume");
	ASSERT_EQ((int)vol.size(), num_agents);
	for (int i = 0; i < num_agents; ++i)
		EXPECT_DOUBLE_EQ(static_cast<double>(vol[i]), (i + 1) * 100.0);

	for (index_t s = 0; s < m->substrates_count; ++s)
	{
		const std::string sname = m->substrates_names[s];

		auto sec = read_point_array(vtu_file, (sname + "_secretion_rate").c_str());
		ASSERT_EQ((int)sec.size(), num_agents);
		auto sat = read_point_array(vtu_file, (sname + "_saturation_density").c_str());
		ASSERT_EQ((int)sat.size(), num_agents);
		auto upt = read_point_array(vtu_file, (sname + "_uptake_rate").c_str());
		ASSERT_EQ((int)upt.size(), num_agents);
		auto net = read_point_array(vtu_file, (sname + "_net_export_rate").c_str());
		ASSERT_EQ((int)net.size(), num_agents);
		auto intern = read_point_array(vtu_file, (sname + "_internalized_substrate").c_str());
		ASSERT_EQ((int)intern.size(), num_agents);
		auto frd = read_point_array(vtu_file, (sname + "_fraction_released_at_death").c_str());
		ASSERT_EQ((int)frd.size(), num_agents);
		auto fti = read_point_array(vtu_file, (sname + "_fraction_transferred_when_ingested").c_str());
		ASSERT_EQ((int)fti.size(), num_agents);

		for (int i = 0; i < num_agents; ++i)
		{
			EXPECT_DOUBLE_EQ(static_cast<double>(sec[i]), i + static_cast<int>(s) + 1.0);
			EXPECT_DOUBLE_EQ(static_cast<double>(sat[i]), (i + 1) * (static_cast<int>(s) + 1) * 10.0);
			EXPECT_DOUBLE_EQ(static_cast<double>(upt[i]), (i + 1) * 0.5);
			EXPECT_DOUBLE_EQ(static_cast<double>(net[i]), (i + 1) * 0.1);
			EXPECT_DOUBLE_EQ(static_cast<double>(intern[i]), i * 5.0);
			EXPECT_DOUBLE_EQ(static_cast<double>(frd[i]), 0.5 + i * 0.1);
			EXPECT_DOUBLE_EQ(static_cast<double>(fti[i]), 0.3 + i * 0.05);
		}
	}
}

TEST_F(VtkAgentsSerializerTest, SerializeMultipleTimes)
{
	auto m = create_test_microenvironment();

	// Create an agent
	m->agents->create();

	vtk_agents_serializer serializer(test_output_dir.string(), *m);

	// Serialize multiple times
	for (int i = 0; i < 3; ++i)
	{
		EXPECT_NO_THROW(serializer.serialize(*m, i * 0.1));
	}

	auto vtk_dir = test_output_dir / "vtk_agents";

	// Check that multiple VTK files are created
	EXPECT_TRUE(std::filesystem::exists(vtk_dir / "agents_000000.vtu"));
	EXPECT_TRUE(std::filesystem::exists(vtk_dir / "agents_000001.vtu"));
	EXPECT_TRUE(std::filesystem::exists(vtk_dir / "agents_000002.vtu"));
}

TEST_F(VtkAgentsSerializerTest, PvdFileContainsCorrectEntries)
{
	auto m = create_test_microenvironment();
	m->agents->create();

	vtk_agents_serializer serializer(test_output_dir.string(), *m);

	// Serialize twice
	serializer.serialize(*m, 0.1);
	serializer.serialize(*m, 0.2);

	// Read PVD file content
	auto pvd_file = test_output_dir / "agents.pvd";
	ASSERT_TRUE(std::filesystem::exists(pvd_file));

	// Check that .vtu files exist
	auto vtk_dir = test_output_dir / "vtk_agents";
	EXPECT_TRUE(std::filesystem::exists(vtk_dir / "agents_000000.vtu"));
	EXPECT_TRUE(std::filesystem::exists(vtk_dir / "agents_000001.vtu"));

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
	EXPECT_TRUE(content.find("agents_000000.vtu") != std::string::npos);
	EXPECT_TRUE(content.find("agents_000001.vtu") != std::string::npos);
}

TEST_F(VtkAgentsSerializerTest, AllSubstrateArraysPresent)
{
	auto m = create_test_microenvironment();
	m->agents->create();

	vtk_agents_serializer serializer(test_output_dir.string(), *m);
	serializer.serialize(*m, 0.0);

	// Read VTK file via pugixml
	auto vtk_dir = test_output_dir / "vtk_agents";
	auto vtu_file = vtk_dir / "agents_000000.vtu";

	pugi::xml_document doc;
	ASSERT_TRUE(doc.load_file(vtu_file.string().c_str()));
	auto point_data = doc.child("VTKFile").child("UnstructuredGrid").child("Piece").child("PointData");

	// Check that all expected arrays are present for each substrate
	const std::vector<std::string> substrates = { "O2", "Glucose" };
	const std::vector<std::string> array_suffixes = { "_secretion_rate",
													  "_saturation_density",
													  "_uptake_rate",
													  "_net_export_rate",
													  "_internalized_substrate",
													  "_fraction_released_at_death",
													  "_fraction_transferred_when_ingested" };

	for (const auto& substrate : substrates)
	{
		for (const auto& suffix : array_suffixes)
		{
			const std::string array_name = substrate + suffix;
			bool found = false;
			for (auto da : point_data.children("DataArray"))
				if (std::string(da.attribute("Name").value()) == array_name) { found = true; break; }
			EXPECT_TRUE(found) << "Missing array: " << array_name;
		}
	}

	// Also check volume array
	bool has_volume = false;
	for (auto da : point_data.children("DataArray"))
		if (std::string(da.attribute("Name").value()) == "volume") { has_volume = true; break; }
	EXPECT_TRUE(has_volume);
}

TEST_F(VtkAgentsSerializerTest, IntegrationWithMicroenvironment)
{
	auto m = create_test_microenvironment();

	// Create agents
	const int num_agents = 3;
	for (int i = 0; i < num_agents; ++i)
	{
		m->agents->create();
	}

	auto container = std::dynamic_pointer_cast<agent_container>(m->agents);
	auto& base_data = std::get<0>(container->agent_datas);
	auto& biofvm_data = std::get<1>(container->agent_datas);

	// Set some agent data
	for (int i = 0; i < num_agents; ++i)
	{
		base_data->positions[i * 3 + 0] = 10.0 + i * 5.0;
		base_data->positions[i * 3 + 1] = 20.0 + i * 5.0;
		base_data->positions[i * 3 + 2] = 30.0 + i * 5.0;
		biofvm_data->volumes[i] = 100.0 * (i + 1);
	}

	// Initialize densities
	m->solver->solve(*m, 1);

	// Serialize through microenvironment's serialize_state method
	EXPECT_NO_THROW(m->serialize_state(0.0));

	// Verify both microenvironment and agent files are created
	auto microenv_pvd = std::filesystem::path("output") / "microenvironment.pvd";
	auto agents_pvd = std::filesystem::path("output") / "agents.pvd";

	EXPECT_TRUE(std::filesystem::exists(microenv_pvd));
	EXPECT_TRUE(std::filesystem::exists(agents_pvd));

	// Verify agent VTU file
	auto vtk_dir = std::filesystem::path("output") / "vtk_agents";
	auto vtu_file = vtk_dir / "agents_000000.vtu";
	EXPECT_TRUE(std::filesystem::exists(vtu_file));

	// Read and verify the file via pugixml
	pugi::xml_document doc;
	ASSERT_TRUE(doc.load_file(vtu_file.string().c_str()));
	auto piece = doc.child("VTKFile").child("UnstructuredGrid").child("Piece");
	EXPECT_EQ(std::stoi(piece.attribute("NumberOfPoints").value()), num_agents);
}
