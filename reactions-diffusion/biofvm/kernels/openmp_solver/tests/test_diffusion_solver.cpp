#include <biofvm/microenvironment.h>
#include <gtest/gtest.h>

#include "diffusion_solver.h"

using namespace physicore;
using namespace physicore::biofvm;

using namespace physicore::biofvm::kernels::openmp_solver;

namespace {
std::unique_ptr<microenvironment> default_microenv(cartesian_mesh mesh)
{
	const real_t timestep = 5;
	const index_t substrates_count = 2;

	auto diff_coefs = std::make_unique<real_t[]>(2);
	diff_coefs[0] = 4;
	diff_coefs[1] = 2;
	auto decay_rates = std::make_unique<real_t[]>(2);
	decay_rates[0] = 5;
	decay_rates[1] = 3;

	auto initial_conds = std::make_unique<real_t[]>(2);
	initial_conds[0] = 1;
	initial_conds[1] = 1;

	auto m = std::make_unique<microenvironment>(mesh, substrates_count, timestep);
	m->diffusion_coefficients = std::move(diff_coefs);
	m->decay_rates = std::move(decay_rates);
	m->initial_conditions = std::move(initial_conds);

	return m;
}

std::unique_ptr<microenvironment> biorobots_microenv(cartesian_mesh mesh)
{
	const real_t timestep = 0.01;
	const index_t substrates_count = 2;

	auto diff_coefs = std::make_unique<real_t[]>(2);
	diff_coefs[0] = 1000;
	diff_coefs[1] = 1000;
	auto decay_rates = std::make_unique<real_t[]>(2);
	decay_rates[0] = 0.1;
	decay_rates[1] = 0.4;

	auto initial_conds = std::make_unique<real_t[]>(2);
	initial_conds[0] = 0;
	initial_conds[1] = 0;

	auto m = std::make_unique<microenvironment>(mesh, substrates_count, timestep);
	m->diffusion_coefficients = std::move(diff_coefs);
	m->decay_rates = std::move(decay_rates);
	m->initial_conditions = std::move(initial_conds);

	return m;
}
} // namespace

TEST(DiffusionSolverTest, Uniform1D)
{
	const cartesian_mesh mesh(1, { 0, 0, 0 }, { 80, 0, 0 }, { 20, 0, 0 });

	auto m = default_microenv(mesh);

	diffusion_solver solver;

	solver.prepare(*m, 1);
	solver.initialize();

#pragma omp parallel
	solver.solve();

	auto dens_l = solver.get_substrates_layout<1>();
	real_t* densities = solver.get_substrates_pointer();

	noarr::traverser(dens_l).for_dims<'x'>([&](auto s) {
		auto l = dens_l ^ noarr::fix(s);

		EXPECT_NEAR(l | noarr::get_at<'s'>(densities, 0), 0.03846154, 1e-6);
		EXPECT_NEAR(l | noarr::get_at<'s'>(densities, 1), 0.0625, 1e-6);
	});
}

TEST(DiffusionSolverTest, Uniform2D)
{
	const cartesian_mesh mesh(2, { 0, 0, 0 }, { 800, 800, 0 }, { 20, 20, 0 });

	auto m = default_microenv(mesh);

	diffusion_solver solver;

	solver.prepare(*m, 1);
	solver.initialize();

#pragma omp parallel
	solver.solve();

	auto dens_l = solver.get_substrates_layout<2>();
	real_t* densities = solver.get_substrates_pointer();

	noarr::traverser(dens_l).for_dims<'x', 'y'>([&](auto s) {
		auto l = dens_l ^ noarr::fix(s);

		EXPECT_NEAR(l | noarr::get_at<'s'>(densities, 0), 0.0054869675, 1e-6);
		EXPECT_NEAR(l | noarr::get_at<'s'>(densities, 1), 0.013840831, 1e-6);
	});
}

TEST(DiffusionSolverTest, Uniform3D)
{
	const cartesian_mesh mesh(3, { 0, 0, 0 }, { 200, 200, 200 }, { 20, 20, 20 });

	auto m = default_microenv(mesh);

	diffusion_solver solver;

	solver.prepare(*m, 1);
	solver.initialize();

#pragma omp parallel
	solver.solve();

	auto dens_l = solver.get_substrates_layout();
	real_t* densities = solver.get_substrates_pointer();

	noarr::traverser(dens_l).for_dims<'x', 'y', 'z'>([&](auto s) {
		auto l = dens_l ^ noarr::fix(s);

		EXPECT_NEAR(l | noarr::get_at<'s'>(densities, 0), 0.0012299563, 1e-6);
		EXPECT_NEAR(l | noarr::get_at<'s'>(densities, 1), 0.0046296306, 1e-6);
	});
}

TEST(DiffusionSolverTest, Random1D)
{
	cartesian_mesh mesh(1, { 0, 0, 0 }, { 60, 0, 0 }, { 20, 20, 20 });

	auto m = biorobots_microenv(mesh);

	diffusion_solver solver;

	solver.prepare(*m, 1);
	solver.initialize();

	auto dens_l = solver.get_substrates_layout<1>();
	real_t* densities = solver.get_substrates_pointer();

	// fill with random values
	for (index_t s = 0; s < m->substrates_count; ++s)
		for (index_t x = 0; x < mesh.grid_shape[0]; ++x)
		{
			const index_t index = s + x * m->substrates_count;
			(dens_l | noarr::get_at<'x', 's'>(densities, x, s)) = static_cast<real_t>(index);
		}

#pragma omp parallel
	solver.solve();

	std::vector<double> expected = {
		0.0486842592, 1.0444132121, 1.9980019980, 2.9880478088, 3.9473197368, 4.9316824055
	};

	for (index_t s = 0; s < m->substrates_count; ++s)
		for (index_t x = 0; x < mesh.grid_shape[0]; ++x)
		{
			const index_t index = s + x * m->substrates_count;

			EXPECT_NEAR((dens_l | noarr::get_at<'x', 's'>(densities, x, s)), expected[index], 1e-6);
		}
}

TEST(DiffusionSolverTest, Random2D)
{
	cartesian_mesh mesh(2, { 0, 0, 0 }, { 60, 60, 0 }, { 20, 20, 20 });

	auto m = biorobots_microenv(mesh);

	diffusion_solver solver;

	solver.prepare(*m, 1);
	solver.initialize();

	auto dens_l = solver.get_substrates_layout<2>();
	real_t* densities = solver.get_substrates_pointer();

	// fill with random values
	for (index_t s = 0; s < m->substrates_count; ++s)
		for (index_t x = 0; x < mesh.grid_shape[0]; ++x)
			for (index_t y = 0; y < mesh.grid_shape[1]; ++y)
			{
				const index_t index = s + x * m->substrates_count + y * m->substrates_count * mesh.grid_shape[0];
				(dens_l | noarr::get_at<'x', 'y', 's'>(densities, x, y, s)) = static_cast<real_t>(index);
			}

#pragma omp parallel
	solver.solve();

	std::vector<double> expected = { 0.1948319355,	1.1899772978,  2.1441254507,  3.1335099015,	 4.0934189658,
									 5.0770425053,	6.0427124809,  7.0205751090,  7.9920058,	 8.9641077127,
									 9.9412995111,	10.9076403164, 11.8905930262, 12.8511729202, 13.8398865413,
									 14.7947055239, 15.7891800565, 16.7382381276 };

	for (index_t s = 0; s < m->substrates_count; ++s)
		for (index_t x = 0; x < mesh.grid_shape[0]; ++x)
			for (index_t y = 0; y < mesh.grid_shape[1]; ++y)
			{
				const index_t index = s + x * m->substrates_count + y * m->substrates_count * mesh.grid_shape[0];

				EXPECT_NEAR((dens_l | noarr::get_at<'x', 'y', 's'>(densities, x, y, s)), expected[index], 1e-6);
			}
}

TEST(DiffusionSolverTest, Random3D)
{
	cartesian_mesh mesh(3, { 0, 0, 0 }, { 60, 60, 60 }, { 20, 20, 20 });

	auto m = biorobots_microenv(mesh);

	diffusion_solver solver;

	solver.prepare(*m, 1);
	solver.initialize();

	auto dens_l = solver.get_substrates_layout<3>();
	real_t* densities = solver.get_substrates_pointer();

	// fill with random values
	for (index_t s = 0; s < m->substrates_count; ++s)
		for (index_t x = 0; x < mesh.grid_shape[0]; ++x)
			for (index_t y = 0; y < mesh.grid_shape[1]; ++y)
				for (index_t z = 0; z < mesh.grid_shape[2]; ++z)
				{
					const index_t index = s + x * m->substrates_count + y * m->substrates_count * mesh.grid_shape[0]
										  + z * m->substrates_count * mesh.grid_shape[0] * mesh.grid_shape[1];
					(dens_l | noarr::get_at<'x', 'y', 'z', 's'>(densities, x, y, z, s)) = static_cast<real_t>(index);
				}

#pragma omp parallel
	solver.solve();

	std::vector<double> expected = {
		0.6333066643,  1.6268066007,  2.5825920996,	 3.5703051208,	4.5318775349,  5.5138036408,  6.4811629703,
		7.4573021609,  8.4304484056,  9.4008006809,	 10.3797338410, 11.3442992010, 12.3290192763, 13.2877977210,
		14.2783047117, 15.2312962410, 16.2275901470, 17.1747947611, 18.1768755823, 19.1182932811, 20.1261610177,
		21.0617918012, 22.0754464530, 23.0052903212, 24.0247318884, 24.9487888412, 25.9740173237, 26.8922873613,
		27.9233027591, 28.8357858813, 29.8725881944, 30.7792844014, 31.8218736297, 32.7227829214, 33.7711590651,
		34.6662814414, 35.7204445004, 36.6097799615, 37.6697299358, 38.5532784815, 39.6190153711, 40.4967770016,
		41.5683008064, 42.4402755216, 43.5175862418, 44.3837740416, 45.4668716771, 46.3272725617, 47.4161571125,
		48.2707710817, 49.3654425478, 50.2142696018, 51.3147279832, 52.1577681218
	};

	for (index_t s = 0; s < m->substrates_count; ++s)
		for (index_t x = 0; x < mesh.grid_shape[0]; ++x)
			for (index_t y = 0; y < mesh.grid_shape[1]; ++y)
				for (index_t z = 0; z < mesh.grid_shape[2]; ++z)
				{
					const index_t index = s + x * m->substrates_count + y * m->substrates_count * mesh.grid_shape[0]
										  + z * m->substrates_count * mesh.grid_shape[0] * mesh.grid_shape[1];

					EXPECT_NEAR((dens_l | noarr::get_at<'x', 'y', 'z', 's'>(densities, x, y, z, s)), expected[index],
								1e-6);
				}
}
