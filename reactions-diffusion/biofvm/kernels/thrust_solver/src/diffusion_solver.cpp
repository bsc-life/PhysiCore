#include "diffusion_solver.h"

#include <noarr/structures/extra/traverser.hpp>
#include <noarr/structures/interop/bag.hpp>
#include <noarr/structures_extended.hpp>
#include <thrust/execution_policy.h>

#include "types.h"

using namespace physicore;
using namespace physicore::biofvm::kernels::thrust_solver;

index_t lcm(index_t a, index_t b)
{
	index_t high = std::max(a, b);
	index_t low = std::min(a, b);

	index_t ret = high;

	while (ret % low != 0)
		ret += high;

	return ret;
}


void diffusion_solver::initialize(microenvironment& m)
{
	// Here we try to find the least common multiple of substrates size in bits and size of a vector register (we assume
	// 512 bits). Thanks to this factor, we can reorganize loops in diffusion so they are automatically vectorized.
	// Finding least common multiple is the most optimal, as the loop has no remainder wrt vector registers.
	// But we want to limit it - when substrates size is much higher than 512 bits, multiplying it will not benefit much

	index_t register_size = 512;
	index_t substrates_size = m.substrates_count * sizeof(real_t) * 8;
	index_t multiple = lcm(substrates_size, register_size);

	while (multiple > 10 * register_size && multiple > substrates_size)
		multiple -= substrates_size;

	initialize(m, multiple / substrates_size);
}

void diffusion_solver::initialize(microenvironment& m, index_t substrate_factor)
{
	substrate_factor_ = substrate_factor;

	if (m.mesh.dims >= 1)
		precompute_values(bx_, cx_, ex_, m.mesh.voxel_shape[0], m.mesh.dims, m.mesh.grid_shape[0], m, 1);
	if (m.mesh.dims >= 2)
		precompute_values(by_, cy_, ey_, m.mesh.voxel_shape[1], m.mesh.dims, m.mesh.grid_shape[1], m,
						  substrate_factor_);
	if (m.mesh.dims >= 3)
		precompute_values(bz_, cz_, ez_, m.mesh.voxel_shape[2], m.mesh.dims, m.mesh.grid_shape[2], m,
						  substrate_factor_);

	ns_ = m.substrates_count;
	nx_ = m.mesh.grid_shape[0];
	ny_ = m.mesh.grid_shape[1];
	nz_ = m.mesh.grid_shape[2];

	substrate_densities_.resize((std::size_t)ns_ * nx_ * ny_ * nz_);
	substrate_densities_.shrink_to_fit();

	initial_conditions_.resize(ns_);
	initial_conditions_.shrink_to_fit();
	thrust::copy(m.initial_conditions.get(), m.initial_conditions.get() + ns_, initial_conditions_.begin());

	thrust::for_each(
		thrust::make_counting_iterator<std::size_t>(0),
		thrust::make_counting_iterator((std::size_t)ns_ * nx_ * ny_ * nz_),
		[ns = ns_, densities = substrate_densities_.data().get(),
		 initial_conditions = initial_conditions_.data().get()] PHYSICORE_THRUST_DEVICE_FN(std::size_t voxel_idx) {
			const index_t s = voxel_idx % ns;
			densities[voxel_idx] = initial_conditions[s];
		});
}

void diffusion_solver::precompute_values(device_vector<real_t>& b, device_vector<real_t>& c, device_vector<real_t>& e,
										 index_t shape, index_t dims, index_t n, const microenvironment& m,
										 index_t copies)
{
	if (n == 1) // special case
	{
		b.resize(m.substrates_count * copies);
		b.shrink_to_fit();

		for (index_t x = 0; x < copies; x++)
			for (index_t s = 0; s < m.substrates_count; s++)
				b[x * m.substrates_count + s] = 1 / (1 + m.decay_rates[s] * m.diffusion_timestep / dims);

		return;
	}

	b.resize(n * m.substrates_count * copies);
	e.resize((n - 1) * m.substrates_count * copies);
	c.resize(m.substrates_count * copies);
	b.shrink_to_fit();
	e.shrink_to_fit();
	c.shrink_to_fit();

	auto layout = noarr::scalar<real_t>() ^ noarr::vector<'s'>() ^ noarr::vector<'x'>() ^ noarr::vector<'i'>()
				  ^ noarr::set_length<'i'>(n) ^ noarr::set_length<'x'>(copies)
				  ^ noarr::set_length<'s'>(m.substrates_count);

	auto b_diag = noarr::make_bag(layout, b.data().get());
	auto e_diag = noarr::make_bag(layout, e.data().get());

	// compute c_i
	for (index_t x = 0; x < copies; x++)
		for (index_t s = 0; s < m.substrates_count; s++)
			c[x * m.substrates_count + s] = -m.diffusion_timestep * m.diffusion_coefficients[s] / (shape * shape);

	// compute b_i
	{
		std::array<index_t, 2> indices = { 0, n - 1 };

		for (index_t i : indices)
			for (index_t x = 0; x < copies; x++)
				for (index_t s = 0; s < m.substrates_count; s++)
					b_diag.at<'i', 'x', 's'>(i, x, s) =
						1 + m.decay_rates[s] * m.diffusion_timestep / dims
						+ m.diffusion_timestep * m.diffusion_coefficients[s] / (shape * shape);

		for (index_t i = 1; i < n - 1; i++)
			for (index_t x = 0; x < copies; x++)
				for (index_t s = 0; s < m.substrates_count; s++)
					b_diag.at<'i', 'x', 's'>(i, x, s) =
						1 + m.decay_rates[s] * m.diffusion_timestep / dims
						+ 2 * m.diffusion_timestep * m.diffusion_coefficients[s] / (shape * shape);
	}

	// compute b_i' and e_i
	{
		for (index_t x = 0; x < copies; x++)
			for (index_t s = 0; s < m.substrates_count; s++)
				b_diag.at<'i', 'x', 's'>(0, x, s) = 1 / b_diag.at<'i', 'x', 's'>(0, x, s);

		for (index_t i = 1; i < n; i++)
			for (index_t x = 0; x < copies; x++)
				for (index_t s = 0; s < m.substrates_count; s++)
				{
					b_diag.at<'i', 'x', 's'>(i, x, s) =
						1
						/ (b_diag.at<'i', 'x', 's'>(i, x, s)
						   - c[x * m.substrates_count + s] * c[x * m.substrates_count + s]
								 * b_diag.at<'i', 'x', 's'>(i - 1, x, s));

					e_diag.at<'i', 'x', 's'>(i - 1, x, s) =
						c[x * m.substrates_count + s] * b_diag.at<'i', 'x', 's'>(i - 1, x, s);
				}
	}
}

template <char swipe_dim, typename density_layout_t, typename index_t>
static constexpr void solve_slice(real_t* _CCCL_RESTRICT densities, const real_t* _CCCL_RESTRICT b,
								  const real_t* _CCCL_RESTRICT c, const real_t* _CCCL_RESTRICT e,
								  const density_layout_t dens_l, index_t s)
{
	const index_t substrates_count = dens_l | noarr::get_length<'s'>();
	const index_t n = dens_l | noarr::get_length<swipe_dim>();

	auto diag_l = noarr::scalar<real_t>() ^ noarr::vector<'s'>(substrates_count) ^ noarr::vector<'i'>(n);

	for (index_t i = 1; i < n; i++)
	{
		(dens_l | noarr::get_at<swipe_dim, 's'>(densities, i, s)) =
			(dens_l | noarr::get_at<swipe_dim, 's'>(densities, i, s))
			- (diag_l | noarr::get_at<'i', 's'>(e, i - 1, s))
				  * (dens_l | noarr::get_at<swipe_dim, 's'>(densities, i - 1, s));
	}

	(dens_l | noarr::get_at<swipe_dim, 's'>(densities, n - 1, s)) =
		(dens_l | noarr::get_at<swipe_dim, 's'>(densities, n - 1, s)) * (diag_l | noarr::get_at<'i', 's'>(b, n - 1, s));

	for (index_t i = n - 2; i >= 0; i--)
	{
		(dens_l | noarr::get_at<swipe_dim, 's'>(densities, i, s)) =
			((dens_l | noarr::get_at<swipe_dim, 's'>(densities, i, s))
			 - c[s] * (dens_l | noarr::get_at<swipe_dim, 's'>(densities, i + 1, s)))
			* (diag_l | noarr::get_at<'i', 's'>(b, i, s));
	}
}

real_t* diffusion_solver::get_substrates_pointer() { return substrate_densities_.data().get(); }

void diffusion_solver::solve()
{
	auto dens_l = get_substrates_layout();

	std::size_t x_work = (std::size_t)ns_ * ny_ * nz_;

	// swipe x
	thrust::for_each(
		thrust::device, thrust::make_counting_iterator<std::size_t>(0), thrust::make_counting_iterator(x_work),
		[dens_l, densities = substrate_densities_.data().get(), b = bx_.data().get(), c = cx_.data().get(),
		 e = ex_.data().get()] PHYSICORE_THRUST_DEVICE_FN(std::size_t voxel_idx) {
			const index_t s_len = dens_l | noarr::get_length<'s'>();
			const index_t y_len = dens_l | noarr::get_length<'y'>();

			const index_t s = voxel_idx % s_len;
			voxel_idx /= s_len;
			const index_t y = voxel_idx % y_len;
			const index_t z = voxel_idx / y_len;

			solve_slice<'x'>(densities, b, c, e, dens_l ^ noarr::fix<'y'>(y) ^ noarr::fix<'z'>(z), (sindex_t)s);
		});

	if (ny_ != 1)
	{
		std::size_t y_work = (std::size_t)ns_ * nx_ * nz_;

		// swipe y
		thrust::for_each(
			thrust::device, thrust::make_counting_iterator<std::size_t>(0), thrust::make_counting_iterator(y_work),
			[dens_l, densities = substrate_densities_.data().get(), b = by_.data().get(), c = cy_.data().get(),
			 e = ey_.data().get()] PHYSICORE_THRUST_DEVICE_FN(std::size_t voxel_idx) {
				const index_t s_len = dens_l | noarr::get_length<'s'>();
				const index_t x_len = dens_l | noarr::get_length<'x'>();

				const index_t s = voxel_idx % s_len;
				voxel_idx /= s_len;
				const index_t x = voxel_idx % x_len;
				const index_t z = voxel_idx / x_len;

				solve_slice<'y'>(densities, b, c, e, dens_l ^ noarr::fix<'x'>(x) ^ noarr::fix<'z'>(z), (sindex_t)s);
			});
	}

	if (nz_ != 1)
	{
		std::size_t z_work = (std::size_t)ns_ * nx_ * ny_;

		// swipe y
		thrust::for_each(
			thrust::device, thrust::make_counting_iterator<std::size_t>(0), thrust::make_counting_iterator(z_work),
			[dens_l, densities = substrate_densities_.data().get(), b = bz_.data().get(), c = cz_.data().get(),
			 e = ez_.data().get()] PHYSICORE_THRUST_DEVICE_FN(std::size_t voxel_idx) {
				const index_t s_len = dens_l | noarr::get_length<'s'>();
				const index_t x_len = dens_l | noarr::get_length<'x'>();

				const index_t s = voxel_idx % s_len;
				voxel_idx /= s_len;
				const index_t x = voxel_idx % x_len;
				const index_t y = voxel_idx / x_len;

				solve_slice<'z'>(densities, b, c, e, dens_l ^ noarr::fix<'x'>(x) ^ noarr::fix<'y'>(y), (sindex_t)s);
			});
	}
}
