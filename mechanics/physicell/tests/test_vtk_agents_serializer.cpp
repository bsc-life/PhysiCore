#include <vtkPointData.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <gmock/gmock.h>

#include "vtk_agents_serializer.h"

using namespace physicore;
using namespace physicore::mechanics::physicell;

namespace {

std::unique_ptr<environment> make_environment(index_t dims, index_t agent_types, index_t substrates)
{
	auto base = std::make_unique<base_agent_data>(dims);
	auto data = std::make_unique<mechanical_agent_data>(*base, agent_types, substrates);
	return std::make_unique<environment>(
		cartesian_mesh { dims, { -100, -100, -100 }, { 100, 100, 100 }, { 20, 20, 20 } }, agent_types, substrates, 0.1);
}

std::unique_ptr<environment> make_single_agent_environment()
{
	auto env = make_environment(3, 2, 2);
	env->agents->create();

	auto* mech_agents = dynamic_cast<mechanical_agent_container*>(env->agents.get());

	auto& base_data = *std::get<0>(mech_agents->agent_datas);
	auto& mech_data = *std::get<1>(mech_agents->agent_datas);

	base_data.positions = { 1.0, 2.0, 3.0 };

	mech_data.velocity = { 0.1, 0.2, 0.3 };
	mech_data.previous_velocity = { 0.4, 0.5, 0.6 };
	mech_data.radius[0] = 4.5;

	mech_data.mechanics_data.cell_cell_adhesion_strength[0] = 0.2;
	mech_data.mechanics_data.cell_BM_adhesion_strength[0] = 0.3;
	mech_data.mechanics_data.cell_cell_repulsion_strength[0] = 1.1;
	mech_data.mechanics_data.cell_BM_repulsion_strength[0] = 0.9;
	mech_data.mechanics_data.relative_maximum_adhesion_distance[0] = 1.25;
	mech_data.mechanics_data.maximum_number_of_attachments[0] = 5;
	mech_data.mechanics_data.attachment_elastic_constant[0] = 0.6;
	mech_data.mechanics_data.attachment_rate[0] = 0.7;
	mech_data.mechanics_data.detachment_rate[0] = 0.8;
	mech_data.mechanics_data.cell_adhesion_affinities = { 0.25, 0.5 };

	mech_data.motility_data.is_motile[0] = 1;
	mech_data.motility_data.persistence_time[0] = 5.5;
	mech_data.motility_data.migration_speed[0] = 6.6;
	mech_data.motility_data.migration_bias_direction = { 0.0, 1.0, 0.0 };
	mech_data.motility_data.migration_bias[0] = 0.7;
	mech_data.motility_data.motility_vector = { 1.0, 2.0, 3.0 };
	mech_data.motility_data.restrict_to_2d[0] = 0;
	mech_data.motility_data.chemotaxis_index[0] = 1;
	mech_data.motility_data.chemotaxis_direction[0] = 2;
	mech_data.motility_data.chemotactic_sensitivities = { 0.01, 0.02 };

	mech_data.state_data.simple_pressure[0] = 12.0;
	mech_data.state_data.agent_type_index[0] = 1;
	mech_data.state_data.is_movable[0] = 1;

	return env;
}

vtkRealArray* arr(vtkPointData* point_data, const char* name)
{
	auto* array = vtkRealArray::SafeDownCast(point_data->GetArray(name));
	if (!array)
	{
		ADD_FAILURE() << "Missing array: " << name;
	}
	return array;
}

void expect_scalar(vtkPointData* point_data, const char* name, double value)
{
	auto* array = arr(point_data, name);
	ASSERT_NE(array, nullptr);
	EXPECT_DOUBLE_EQ(array->GetTuple1(0), value);
}

void expect_vec3(vtkPointData* point_data, const char* name, const std::array<double, 3>& expected)
{
	std::array<double, 3> tuple {};
	auto* array = arr(point_data, name);
	ASSERT_NE(array, nullptr);
	array->GetTuple(0, tuple.data());
	EXPECT_DOUBLE_EQ(tuple[0], expected[0]);
	EXPECT_DOUBLE_EQ(tuple[1], expected[1]);
	EXPECT_DOUBLE_EQ(tuple[2], expected[2]);
}

void expect_point_at(vtkUnstructuredGrid* grid, const std::array<double, 3>& expected)
{
	std::array<double, 3> pos {};
	grid->GetPoint(0, pos.data());
	EXPECT_DOUBLE_EQ(pos[0], expected[0]);
	EXPECT_DOUBLE_EQ(pos[1], expected[1]);
	EXPECT_DOUBLE_EQ(pos[2], expected[2]);
}

void assert_mechanics(vtkPointData* point_data)
{
	expect_scalar(point_data, "radius", 4.5);
	expect_scalar(point_data, "cell_cell_adhesion_strength", 0.2);
	expect_scalar(point_data, "cell_BM_adhesion_strength", 0.3);
	expect_scalar(point_data, "cell_cell_repulsion_strength", 1.1);
	expect_scalar(point_data, "cell_BM_repulsion_strength", 0.9);
	expect_scalar(point_data, "relative_maximum_adhesion_distance", 1.25);
	expect_scalar(point_data, "maximum_number_of_attachments", 5.0);
	expect_scalar(point_data, "attachment_elastic_constant", 0.6);
	expect_scalar(point_data, "attachment_rate", 0.7);
	expect_scalar(point_data, "detachment_rate", 0.8);
}

void assert_motility(vtkPointData* point_data)
{
	expect_scalar(point_data, "is_motile", 1.0);
	expect_scalar(point_data, "persistence_time", 5.5);
	expect_scalar(point_data, "migration_speed", 6.6);
	expect_scalar(point_data, "migration_bias", 0.7);
	expect_scalar(point_data, "restrict_to_2d", 0.0);
	expect_scalar(point_data, "chemotaxis_index", 1.0);
	expect_scalar(point_data, "chemotaxis_direction", 2.0);
}

void assert_state(vtkPointData* point_data)
{
	expect_scalar(point_data, "simple_pressure", 12.0);
	expect_scalar(point_data, "cell_definition_index", 1.0);
	expect_scalar(point_data, "is_movable", 1.0);
}

void assert_vectors(vtkPointData* point_data)
{
	expect_vec3(point_data, "velocity", { 0.1, 0.2, 0.3 });
	expect_vec3(point_data, "previous_velocity", { 0.4, 0.5, 0.6 });
	expect_vec3(point_data, "migration_bias_direction", { 0.0, 1.0, 0.0 });
	expect_vec3(point_data, "motility_vector", { 1.0, 2.0, 3.0 });
}

void assert_affinities(vtkPointData* point_data)
{
	expect_scalar(point_data, "immune_cell_adhesion_affinity", 0.25);
	expect_scalar(point_data, "cell_type_1_cell_adhesion_affinity", 0.5);
	expect_scalar(point_data, "O2_chemotactic_sensitivity", 0.01);
	expect_scalar(point_data, "substrate_1_chemotactic_sensitivity", 0.02);
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
	auto e = make_environment(3, 1, 1);

	EXPECT_NO_THROW({ const vtk_agents_serializer serializer(test_output_dir.string(), *e, { "O2", "" }); });

	auto vtk_dir = test_output_dir / "vtk_mechanics_agents";
	EXPECT_TRUE(std::filesystem::exists(vtk_dir));
}

TEST_F(VtkMechanicsAgentsSerializerTest, SerializeWithNoAgentsCreatesFiles)
{
	auto e = make_environment(3, 1, 1);
	vtk_agents_serializer serializer(test_output_dir.string(), *e, { "O2", "" });

	EXPECT_NO_THROW(serializer.serialize(*e, 0.0));

	auto vtk_dir = test_output_dir / "vtk_mechanics_agents";
	EXPECT_TRUE(std::filesystem::exists(vtk_dir / "mechanics_agents_000000.vtu"));
	EXPECT_TRUE(std::filesystem::exists(test_output_dir / "mechanics_agents.pvd"));
}

TEST_F(VtkMechanicsAgentsSerializerTest, SerializeSingleAgentWritesExpectedArrays)
{
	auto e = make_single_agent_environment();

	std::vector<std::string> const substrate_names { "O2", "" };
	e->agent_types.emplace_back();
	e->agent_types.emplace_back();
	e->agent_types[0].name = "immune";
	e->agent_types[1].name = "";

	vtk_agents_serializer serializer(test_output_dir.string(), *e, substrate_names);
	serializer.serialize(*e, 0.0);

	auto vtk_dir = test_output_dir / "vtk_mechanics_agents";
	auto vtu_file = vtk_dir / "mechanics_agents_000000.vtu";

	auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	reader->SetFileName(vtu_file.string().c_str());
	reader->Update();

	auto* grid = reader->GetOutput();
	ASSERT_NE(grid, nullptr);
	EXPECT_EQ(grid->GetNumberOfPoints(), 1);
	expect_point_at(grid, { 1.0, 2.0, 3.0 });

	auto* point_data = grid->GetPointData();
	ASSERT_NE(point_data, nullptr);

	assert_mechanics(point_data);
	assert_motility(point_data);
	assert_state(point_data);
	assert_vectors(point_data);
	assert_affinities(point_data);
}

TEST_F(VtkMechanicsAgentsSerializerTest, SerializeMultipleAgentsWritesAllData)
{
	const int agent_count = 3;
	auto e = make_environment(3, 2, 2);

	for (int i = 0; i < agent_count; ++i)
	{
		auto* agent = e->agents->create();

		agent->position()[0] = 1.0 + i;
		agent->position()[1] = 2.0 + i;
		agent->position()[2] = 3.0 + i;

		agent->radius() = 1.0 + i;
		agent->agent_type_index() = static_cast<index_t>(i);

		agent->velocity()[0] = 0.1 * (i + 1);
		agent->velocity()[1] = 0.2 * (i + 1);
		agent->velocity()[2] = 0.3 * (i + 1);

		agent->chemotactic_sensitivities()[0] = 0.01 * (i + 1);
		agent->chemotactic_sensitivities()[1] = 0.02 * (i + 1);

		agent->cell_adhesion_affinities()[0] = 0.1 * (i + 1);
		agent->cell_adhesion_affinities()[1] = 0.2 * (i + 1);
	}

	std::vector<std::string> const substrate_names { "S1", "S2" };
	e->agent_types.emplace_back();
	e->agent_types.emplace_back();
	e->agent_types[0].name = "typeA";
	e->agent_types[1].name = "typeB";

	vtk_agents_serializer serializer(test_output_dir.string(), *e, substrate_names);
	serializer.serialize(*e, 0.0);

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
		auto* array = vtkRealArray::SafeDownCast(point_data->GetArray(name));
		if (!array)
		{
			ADD_FAILURE() << "Missing array: " << name;
			return nullptr;
		}
		return array;
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
	auto e = make_environment(3, 1, 1);
	auto* agent = e->agents->create();

	agent->position()[0] = 0.0;
	agent->position()[1] = 0.0;
	agent->position()[2] = 0.0;

	vtk_agents_serializer serializer(test_output_dir.string(), *e, {});

	serializer.serialize(*e, 0.0);
	serializer.serialize(*e, 0.1);

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
