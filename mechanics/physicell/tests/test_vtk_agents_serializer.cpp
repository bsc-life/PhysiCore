#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <vector>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <gtest/gtest.h>

#include "physicell/mechanical_agent_container.h"
#include "physicell/vtk_agents_serializer.h"

using namespace physicore;
using namespace physicore::mechanics::physicell;

namespace {

mechanical_agent_container make_container(int dims, int agent_types, int substrates)
{
	auto base = std::make_unique<base_agent_data>(dims);
	auto data = std::make_unique<agent_data>(*base, agent_types, substrates);
	return mechanical_agent_container(std::move(base), std::move(data));
}

} // namespace

class VtkMechanicsAgentsSerializerTest : public ::testing::Test
{
protected:
	std::filesystem::path test_output_dir;

	void SetUp() override
	{
		test_output_dir = std::filesystem::temp_directory_path() / "vtk_mechanics_agents_serializer_test";
		std::filesystem::remove_all(test_output_dir);
		std::filesystem::create_directories(test_output_dir);
	}

	void TearDown() override
	{
		if (std::filesystem::exists(test_output_dir))
		{
			std::filesystem::remove_all(test_output_dir);
		}
	}
};

TEST_F(VtkMechanicsAgentsSerializerTest, ConstructorCreatesOutputDirectories)
{
	auto container = make_container(3, 1, 1);

	EXPECT_NO_THROW({ vtk_agents_serializer serializer(test_output_dir.string(), container); });

	auto vtk_dir = test_output_dir / "vtk_mechanics_agents";
	EXPECT_TRUE(std::filesystem::exists(vtk_dir));
}

TEST_F(VtkMechanicsAgentsSerializerTest, SerializeWithNoAgentsCreatesFiles)
{
	auto container = make_container(3, 1, 1);
	vtk_agents_serializer serializer(test_output_dir.string(), container);

	EXPECT_NO_THROW(serializer.serialize(0.0));

	auto vtk_dir = test_output_dir / "vtk_mechanics_agents";
	EXPECT_TRUE(std::filesystem::exists(vtk_dir / "mechanics_agents_000000.vtu"));
	EXPECT_TRUE(std::filesystem::exists(test_output_dir / "mechanics_agents.pvd"));
}

TEST_F(VtkMechanicsAgentsSerializerTest, SerializeSingleAgentWritesExpectedArrays)
{
	auto container = make_container(3, 2, 2);
	container.create();

	auto& base_data = *std::get<0>(container.agent_datas);
	auto& mech_data = *std::get<1>(container.agent_datas);

	base_data.positions = { 1.0, 2.0, 3.0 };

	mech_data.velocity = { 0.1, 0.2, 0.3 };
	mech_data.previous_velocity = { 0.4, 0.5, 0.6 };
	mech_data.radius[0] = 4.5;

	mech_data.mechanics.cell_cell_adhesion_strength[0] = 0.2;
	mech_data.mechanics.cell_BM_adhesion_strength[0] = 0.3;
	mech_data.mechanics.cell_cell_repulsion_strength[0] = 1.1;
	mech_data.mechanics.cell_BM_repulsion_strength[0] = 0.9;
	mech_data.mechanics.relative_maximum_adhesion_distance[0] = 1.25;
	mech_data.mechanics.maximum_number_of_attachments[0] = 5;
	mech_data.mechanics.attachment_elastic_constant[0] = 0.6;
	mech_data.mechanics.attachment_rate[0] = 0.7;
	mech_data.mechanics.detachment_rate[0] = 0.8;
	mech_data.mechanics.cell_adhesion_affinities = { 0.25, 0.5 };

	mech_data.motility.is_motile[0] = 1;
	mech_data.motility.persistence_time[0] = 5.5;
	mech_data.motility.migration_speed[0] = 6.6;
	mech_data.motility.migration_bias_direction = { 0.0, 1.0, 0.0 };
	mech_data.motility.migration_bias[0] = 0.7;
	mech_data.motility.motility_vector = { 1.0, 2.0, 3.0 };
	mech_data.motility.restrict_to_2d[0] = 0;
	mech_data.motility.chemotaxis_index[0] = 1;
	mech_data.motility.chemotaxis_direction[0] = 2;
	mech_data.motility.chemotactic_sensitivities = { 0.01, 0.02 };

	mech_data.state.orientation = { 0.1, 0.2, 0.3 };
	mech_data.state.simple_pressure[0] = 12.0;
	mech_data.state.agent_type_index[0] = 1;
	mech_data.state.is_movable[0] = 1;

	vtk_agents_serializer serializer(test_output_dir.string(), container, { "O2", "" }, { "immune", "" });
	serializer.serialize(0.0);

	auto vtk_dir = test_output_dir / "vtk_mechanics_agents";
	auto vtu_file = vtk_dir / "mechanics_agents_000000.vtu";

	auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	reader->SetFileName(vtu_file.string().c_str());
	reader->Update();

	auto* grid = reader->GetOutput();
	ASSERT_NE(grid, nullptr);
	EXPECT_EQ(grid->GetNumberOfPoints(), 1);

	std::array<double, 3> pos {};
	grid->GetPoint(0, pos.data());
	EXPECT_DOUBLE_EQ(pos[0], 1.0);
	EXPECT_DOUBLE_EQ(pos[1], 2.0);
	EXPECT_DOUBLE_EQ(pos[2], 3.0);

	auto* point_data = grid->GetPointData();
	ASSERT_NE(point_data, nullptr);

	auto get_arr = [&](const char* name) -> vtkRealArray* {
		auto* arr = point_data->GetArray(name);
		ASSERT_NE(arr, nullptr) << "Missing array: " << name;
		return arr;
	};

	EXPECT_DOUBLE_EQ(get_arr("radius")->GetTuple1(0), 4.5);

	EXPECT_DOUBLE_EQ(get_arr("cell_cell_adhesion_strength")->GetTuple1(0), 0.2);
	EXPECT_DOUBLE_EQ(get_arr("cell_BM_adhesion_strength")->GetTuple1(0), 0.3);
	EXPECT_DOUBLE_EQ(get_arr("cell_cell_repulsion_strength")->GetTuple1(0), 1.1);
	EXPECT_DOUBLE_EQ(get_arr("cell_BM_repulsion_strength")->GetTuple1(0), 0.9);
	EXPECT_DOUBLE_EQ(get_arr("relative_maximum_adhesion_distance")->GetTuple1(0), 1.25);
	EXPECT_DOUBLE_EQ(get_arr("maximum_number_of_attachments")->GetTuple1(0), 5.0);
	EXPECT_DOUBLE_EQ(get_arr("attachment_elastic_constant")->GetTuple1(0), 0.6);
	EXPECT_DOUBLE_EQ(get_arr("attachment_rate")->GetTuple1(0), 0.7);
	EXPECT_DOUBLE_EQ(get_arr("detachment_rate")->GetTuple1(0), 0.8);

	EXPECT_DOUBLE_EQ(get_arr("is_motile")->GetTuple1(0), 1.0);
	EXPECT_DOUBLE_EQ(get_arr("persistence_time")->GetTuple1(0), 5.5);
	EXPECT_DOUBLE_EQ(get_arr("migration_speed")->GetTuple1(0), 6.6);
	EXPECT_DOUBLE_EQ(get_arr("migration_bias")->GetTuple1(0), 0.7);
	EXPECT_DOUBLE_EQ(get_arr("restrict_to_2d")->GetTuple1(0), 0.0);
	EXPECT_DOUBLE_EQ(get_arr("chemotaxis_index")->GetTuple1(0), 1.0);
	EXPECT_DOUBLE_EQ(get_arr("chemotaxis_direction")->GetTuple1(0), 2.0);
	EXPECT_DOUBLE_EQ(get_arr("simple_pressure")->GetTuple1(0), 12.0);
	EXPECT_DOUBLE_EQ(get_arr("cell_definition_index")->GetTuple1(0), 1.0);
	EXPECT_DOUBLE_EQ(get_arr("is_movable")->GetTuple1(0), 1.0);

	std::array<double, 3> tuple {};
	get_arr("velocity")->GetTuple(0, tuple.data());
	EXPECT_DOUBLE_EQ(tuple[0], 0.1);
	EXPECT_DOUBLE_EQ(tuple[1], 0.2);
	EXPECT_DOUBLE_EQ(tuple[2], 0.3);

	get_arr("previous_velocity")->GetTuple(0, tuple.data());
	EXPECT_DOUBLE_EQ(tuple[0], 0.4);
	EXPECT_DOUBLE_EQ(tuple[1], 0.5);
	EXPECT_DOUBLE_EQ(tuple[2], 0.6);

	get_arr("migration_bias_direction")->GetTuple(0, tuple.data());
	EXPECT_DOUBLE_EQ(tuple[0], 0.0);
	EXPECT_DOUBLE_EQ(tuple[1], 1.0);
	EXPECT_DOUBLE_EQ(tuple[2], 0.0);

	get_arr("motility_vector")->GetTuple(0, tuple.data());
	EXPECT_DOUBLE_EQ(tuple[0], 1.0);
	EXPECT_DOUBLE_EQ(tuple[1], 2.0);
	EXPECT_DOUBLE_EQ(tuple[2], 3.0);

	get_arr("orientation")->GetTuple(0, tuple.data());
	EXPECT_DOUBLE_EQ(tuple[0], 0.1);
	EXPECT_DOUBLE_EQ(tuple[1], 0.2);
	EXPECT_DOUBLE_EQ(tuple[2], 0.3);

	EXPECT_DOUBLE_EQ(get_arr("immune_cell_adhesion_affinity")->GetTuple1(0), 0.25);
	EXPECT_DOUBLE_EQ(get_arr("cell_type_1_cell_adhesion_affinity")->GetTuple1(0), 0.5);
	EXPECT_DOUBLE_EQ(get_arr("O2_chemotactic_sensitivity")->GetTuple1(0), 0.01);
	EXPECT_DOUBLE_EQ(get_arr("substrate_1_chemotactic_sensitivity")->GetTuple1(0), 0.02);
}

TEST_F(VtkMechanicsAgentsSerializerTest, SerializeMultipleAgentsWritesAllData)
{
	const int agent_count = 3;
	auto container = make_container(3, 2, 2);
	for (int i = 0; i < agent_count; ++i)
	{
		container.create();
	}

	auto& base_data = *std::get<0>(container.agent_datas);
	auto& mech_data = *std::get<1>(container.agent_datas);

	for (int i = 0; i < agent_count; ++i)
	{
		base_data.positions[i * 3 + 0] = 1.0 + i;
		base_data.positions[i * 3 + 1] = 2.0 + i;
		base_data.positions[i * 3 + 2] = 3.0 + i;

		mech_data.radius[i] = 1.0 + i;
		mech_data.state.agent_type_index[i] = static_cast<index_t>(i);

		mech_data.velocity[i * 3 + 0] = 0.1 * (i + 1);
		mech_data.velocity[i * 3 + 1] = 0.2 * (i + 1);
		mech_data.velocity[i * 3 + 2] = 0.3 * (i + 1);

		mech_data.motility.chemotactic_sensitivities[i * 2 + 0] = 0.01 * (i + 1);
		mech_data.motility.chemotactic_sensitivities[i * 2 + 1] = 0.02 * (i + 1);

		mech_data.mechanics.cell_adhesion_affinities[i * 2 + 0] = 0.1 * (i + 1);
		mech_data.mechanics.cell_adhesion_affinities[i * 2 + 1] = 0.2 * (i + 1);
	}

	vtk_agents_serializer serializer(test_output_dir.string(), container, { "S1", "S2" }, { "typeA", "typeB" });
	serializer.serialize(0.0);

	auto vtu_file = test_output_dir / "vtk_mechanics_agents" / "mechanics_agents_000000.vtu";

	auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	reader->SetFileName(vtu_file.string().c_str());
	reader->Update();

	auto* grid = reader->GetOutput();
	ASSERT_NE(grid, nullptr);
	EXPECT_EQ(grid->GetNumberOfPoints(), agent_count);

	auto* point_data = grid->GetPointData();
	ASSERT_NE(point_data, nullptr);

	auto get_arr = [&](const char* name) -> vtkRealArray* {
		auto* arr = point_data->GetArray(name);
		ASSERT_NE(arr, nullptr) << "Missing array: " << name;
		return arr;
	};

	for (int i = 0; i < agent_count; ++i)
	{
		std::array<double, 3> pos {};
		grid->GetPoint(i, pos.data());
		EXPECT_DOUBLE_EQ(pos[0], 1.0 + i);
		EXPECT_DOUBLE_EQ(pos[1], 2.0 + i);
		EXPECT_DOUBLE_EQ(pos[2], 3.0 + i);

		EXPECT_DOUBLE_EQ(get_arr("radius")->GetTuple1(i), 1.0 + i);
		EXPECT_DOUBLE_EQ(get_arr("cell_definition_index")->GetTuple1(i), i);
		EXPECT_DOUBLE_EQ(get_arr("typeA_cell_adhesion_affinity")->GetTuple1(i), 0.1 * (i + 1));
		EXPECT_DOUBLE_EQ(get_arr("typeB_cell_adhesion_affinity")->GetTuple1(i), 0.2 * (i + 1));
		EXPECT_DOUBLE_EQ(get_arr("S1_chemotactic_sensitivity")->GetTuple1(i), 0.01 * (i + 1));
		EXPECT_DOUBLE_EQ(get_arr("S2_chemotactic_sensitivity")->GetTuple1(i), 0.02 * (i + 1));

		std::array<double, 3> vel {};
		get_arr("velocity")->GetTuple(i, vel.data());
		EXPECT_DOUBLE_EQ(vel[0], 0.1 * (i + 1));
		EXPECT_DOUBLE_EQ(vel[1], 0.2 * (i + 1));
		EXPECT_DOUBLE_EQ(vel[2], 0.3 * (i + 1));
	}
}

TEST_F(VtkMechanicsAgentsSerializerTest, SerializeMultipleTimesAppendsPvd)
{
	auto container = make_container(3, 1, 1);
	container.create();

	auto& base_data = *std::get<0>(container.agent_datas);
	base_data.positions = { 0.0, 0.0, 0.0 };

	vtk_agents_serializer serializer(test_output_dir.string(), container);

	serializer.serialize(0.0);
	serializer.serialize(0.1);

	auto vtk_dir = test_output_dir / "vtk_mechanics_agents";
	EXPECT_TRUE(std::filesystem::exists(vtk_dir / "mechanics_agents_000000.vtu"));
	EXPECT_TRUE(std::filesystem::exists(vtk_dir / "mechanics_agents_000001.vtu"));

	auto pvd_file = test_output_dir / "mechanics_agents.pvd";
	ASSERT_TRUE(std::filesystem::exists(pvd_file));

	std::ifstream file(pvd_file);
	const std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

	EXPECT_NE(content.find("mechanics_agents_000000.vtu"), std::string::npos);
	EXPECT_NE(content.find("mechanics_agents_000001.vtu"), std::string::npos);
	EXPECT_NE(content.find("timestep=\"0.0"), std::string::npos);
	EXPECT_NE(content.find("timestep=\"0.1"), std::string::npos);
}
