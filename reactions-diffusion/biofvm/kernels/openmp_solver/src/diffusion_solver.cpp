#include "diffusion_solver.h"

#include <common/types.h>
#include <hwy/aligned_allocator.h>
#include <hwy/base.h>
#include <noarr/structures/interop/bag.hpp>

#include "solver_utils.h"

using namespace physicore;
using namespace physicore::biofvm::kernels::openmp_solver;

namespace {
auto get_diagonal_layout(const problem_t& problem, index_t n)
{
	return noarr::scalar<real_t>() ^ noarr::vectors<'s', 'i'>(problem.substrates_count, n);
}

auto get_diagonal_layout_c(const problem_t& problem, index_t n, index_t copies)
{
	return noarr::scalar<real_t>() ^ noarr::vectors<'c', 'i'>(problem.substrates_count * copies, n);
}
} // namespace

real_t* diffusion_solver::get_substrates_pointer() { return substrates_.get(); }

const real_t* diffusion_solver::get_substrates_pointer() const { return substrates_.get(); }

void diffusion_solver::precompute_values(std::unique_ptr<real_t[]>& b, std::unique_ptr<real_t[]>& c,
										 std::unique_ptr<real_t[]>& e, index_t shape, index_t dims, index_t n,
										 index_t copies)
{
	b = std::make_unique<real_t[]>(n * problem.substrates_count * copies);
	e = std::make_unique<real_t[]>((n - 1) * problem.substrates_count * copies);
	c = std::make_unique<real_t[]>(problem.substrates_count * copies);

	auto layout = noarr::scalar<real_t>() ^ noarr::vector<'s'>(problem.substrates_count) ^ noarr::vector<'c'>(copies)
				  ^ noarr::vector<'i'>(n);

	auto b_diag = noarr::make_bag(layout, b.get());
	auto e_diag = noarr::make_bag(layout, e.get());

	// compute c_i
	for (index_t x = 0; x < copies; x++)
		for (index_t s = 0; s < problem.substrates_count; s++)
			c[x * problem.substrates_count + s] =
				-1 * -problem.dt * problem.diffusion_coefficients[s] / (real_t)(shape * shape);

	// compute b_i
	{
		for (index_t i = 0; i < n; i++)
			for (index_t x = 0; x < copies; x++)
				for (index_t s = 0; s < problem.substrates_count; s++)
				{
					b_diag.template at<'i', 'c', 's'>(i, x, s) =
						1 + problem.decay_rates[s] * problem.dt / (real_t)dims
						+ 2 * problem.dt * problem.diffusion_coefficients[s] / (real_t)(shape * shape);

					if (i == 0 || i == n - 1)
						b_diag.template at<'i', 'c', 's'>(i, x, s) -=
							problem.dt * problem.diffusion_coefficients[s] / (real_t)(shape * shape);
				}
	}

	// compute b_i' and e_i
	{
		for (index_t copy_idx = 0; copy_idx < copies; copy_idx++)
			for (index_t s = 0; s < problem.substrates_count; s++)
				b_diag.template at<'i', 'c', 's'>(0, copy_idx, s) =
					1 / b_diag.template at<'i', 'c', 's'>(0, copy_idx, s);

		for (index_t i = 1; i < n; i++)
			for (index_t x = 0; x < copies; x++)
				for (index_t s = 0; s < problem.substrates_count; s++)
				{
					b_diag.template at<'i', 'c', 's'>(i, x, s) =
						1
						/ (b_diag.template at<'i', 'c', 's'>(i, x, s)
						   - c[x * problem.substrates_count + s] * c[x * problem.substrates_count + s]
								 * b_diag.template at<'i', 'c', 's'>(i - 1, x, s));

					e_diag.template at<'i', 'c', 's'>(i - 1, x, s) =
						c[x * problem.substrates_count + s] * b_diag.template at<'i', 'c', 's'>(i - 1, x, s);
				}
	}
}


void diffusion_solver::prepare(const microenvironment& m, index_t iterations)
{
	problem = problem_t::construct(m, iterations);

	auto substrates_layout = get_substrates_layout<3>();

	this->substrates_ = hwy::MakeUniqueAlignedArray<real_t>((substrates_layout | noarr::get_size()) / sizeof(real_t));

	// Initialize substrates
	solver_utils::initialize_substrate_constant(substrates_layout, this->substrates_.get(),
												problem.initial_conditions.data());
}


void diffusion_solver::initialize()
{
	substrate_copies_ = (xs_tile_size_ + problem.substrates_count - 1) / problem.substrates_count;

	if (problem.dims >= 1)
		precompute_values(bx_, cx_, ex_, problem.dx, problem.dims, problem.nx, 1);
	if (problem.dims >= 2)
		precompute_values(by_, cy_, ey_, problem.dy, problem.dims, problem.ny, substrate_copies_);
	if (problem.dims >= 3)
		precompute_values(bz_, cz_, ez_, problem.dz, problem.dims, problem.nz, substrate_copies_);
}

namespace {
template <typename index_t, typename real_t, typename density_layout_t, typename diagonal_layout_t>
void solve_slice_x_1d(real_t* HWY_RESTRICT densities, const real_t* HWY_RESTRICT b, const real_t* HWY_RESTRICT c,
					  const real_t* HWY_RESTRICT e, const density_layout_t dens_l, const diagonal_layout_t diag_l)
{
	const index_t substrates_count = dens_l | noarr::get_length<'s'>();
	const index_t n = dens_l | noarr::get_length<'x'>();

	for (index_t i = 1; i < n; i++)
	{
#pragma omp for schedule(static) nowait
		for (index_t s = 0; s < substrates_count; s++)
		{
			(dens_l | noarr::get_at<'x', 's'>(densities, i, s)) =
				(dens_l | noarr::get_at<'x', 's'>(densities, i, s))
				+ (diag_l | noarr::get_at<'i', 's'>(e, i - 1, s))
					  * (dens_l | noarr::get_at<'x', 's'>(densities, i - 1, s));

			// std::cout << i << ": " << (dens_l | noarr::get_at<'x', 's'>(densities, i, s)) << std::endl;
		}
	}

#pragma omp for schedule(static) nowait
	for (index_t s = 0; s < substrates_count; s++)
	{
		(dens_l | noarr::get_at<'x', 's'>(densities, n - 1, s)) =
			(dens_l | noarr::get_at<'x', 's'>(densities, n - 1, s)) * (diag_l | noarr::get_at<'i', 's'>(b, n - 1, s));

		// std::cout << "n-1: " << (dens_l | noarr::get_at<'x', 's'>(densities, n - 1, s)) << std::endl;
	}

	for (index_t i = n - 2; i >= 0; i--)
	{
#pragma omp for schedule(static) nowait
		for (index_t s = 0; s < substrates_count; s++)
		{
			(dens_l | noarr::get_at<'x', 's'>(densities, i, s)) =
				((dens_l | noarr::get_at<'x', 's'>(densities, i, s))
				 + c[s] * (dens_l | noarr::get_at<'x', 's'>(densities, i + 1, s)))
				* (diag_l | noarr::get_at<'i', 's'>(b, i, s));

			// std::cout << i << ": " << (dens_l | noarr::get_at<'x', 's'>(densities, i, s)) << std::endl;
		}
	}
}

template <typename index_t, typename real_t, typename density_layout_t, typename diagonal_layout_t>
void solve_slice_x_2d_and_3d(real_t* HWY_RESTRICT densities, const real_t* HWY_RESTRICT b, const real_t* HWY_RESTRICT c,
							 const real_t* HWY_RESTRICT e, const density_layout_t dens_l,
							 const diagonal_layout_t diag_l)
{
	const index_t substrates_count = dens_l | noarr::get_length<'s'>();
	const index_t n = dens_l | noarr::get_length<'x'>();
	const index_t m = dens_l | noarr::get_length<'m'>();

#pragma omp for schedule(static) nowait
	for (index_t yz = 0; yz < m; yz++)
	{
		for (index_t i = 1; i < n; i++)
		{
			for (index_t s = 0; s < substrates_count; s++)
			{
				(dens_l | noarr::get_at<'m', 'x', 's'>(densities, yz, i, s)) =
					(dens_l | noarr::get_at<'m', 'x', 's'>(densities, yz, i, s))
					+ (diag_l | noarr::get_at<'i', 's'>(e, i - 1, s))
						  * (dens_l | noarr::get_at<'m', 'x', 's'>(densities, yz, i - 1, s));
			}
		}

		for (index_t s = 0; s < substrates_count; s++)
		{
			(dens_l | noarr::get_at<'m', 'x', 's'>(densities, yz, n - 1, s)) =
				(dens_l | noarr::get_at<'m', 'x', 's'>(densities, yz, n - 1, s))
				* (diag_l | noarr::get_at<'i', 's'>(b, n - 1, s));
		}

		for (index_t i = n - 2; i >= 0; i--)
		{
			for (index_t s = 0; s < substrates_count; s++)
			{
				(dens_l | noarr::get_at<'m', 'x', 's'>(densities, yz, i, s)) =
					((dens_l | noarr::get_at<'m', 'x', 's'>(densities, yz, i, s))
					 + c[s] * (dens_l | noarr::get_at<'m', 'x', 's'>(densities, yz, i + 1, s)))
					* (diag_l | noarr::get_at<'i', 's'>(b, i, s));
			}
		}
	}
}

template <typename index_t, typename real_t, typename density_layout_t, typename diagonal_layout_t>
void solve_slice_y_2d(real_t* HWY_RESTRICT densities, const real_t* HWY_RESTRICT b, const real_t* HWY_RESTRICT c_,
					  const real_t* HWY_RESTRICT e, const density_layout_t dens_l, const diagonal_layout_t diag_l,
					  std::size_t s_copies, std::size_t xs_tile_size)
{
	const index_t substrate_count = dens_l | noarr::get_length<'s'>();
	const index_t n = dens_l | noarr::get_length<'y'>();

	auto blocked_dens_l = dens_l ^ noarr::merge_blocks<'x', 's', 'c'>()
						  ^ noarr::into_blocks_static<'c', 'b', 'x', 's'>(substrate_count * s_copies)
						  ^ noarr::into_blocks_static<'s', 'p', 'S', 's'>(xs_tile_size);

	// body
	{
		auto b_dens_l = blocked_dens_l ^ noarr::fix<'b'>(noarr::lit<0>);
		const index_t x_len = b_dens_l | noarr::get_length<'x'>();

#pragma omp for schedule(static) nowait
		for (index_t x = 0; x < x_len; x++)
		{
			// body
			{
				auto bb_dens_l = b_dens_l ^ noarr::fix<'p'>(noarr::lit<0>);
				const index_t S_len = bb_dens_l | noarr::get_length<'S'>();
				const index_t s_len = bb_dens_l | noarr::get_length<'s'>();

				for (index_t S = 0; S < S_len; S++)
				{
					for (index_t i = 1; i < n; i++)
					{
						for (index_t s = 0; s < s_len; s++)
						{
							(bb_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, x, S, s)) =
								(bb_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, x, S, s))
								+ (diag_l | noarr::get_at<'i', 'c'>(e, i - 1, s))
									  * (bb_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i - 1, x, S, s));
						}
					}
					for (index_t s = 0; s < s_len; s++)
					{
						(bb_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, n - 1, x, S, s)) =
							(bb_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, n - 1, x, S, s))
							* (diag_l | noarr::get_at<'i', 'c'>(b, n - 1, s));
					}
					for (index_t i = n - 2; i >= 0; i--)
					{
						for (index_t s = 0; s < s_len; s++)
						{
							(bb_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, x, S, s)) =
								((bb_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, x, S, s))
								 + c_[s] * (bb_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i + 1, x, S, s)))
								* (diag_l | noarr::get_at<'i', 'c'>(b, i, s));
						}
					}
				}
			}

			// remainder
			{
				auto br_dens_l = b_dens_l ^ noarr::fix<'p'>(noarr::lit<1>);
				const index_t s_len = br_dens_l | noarr::get_length<'s'>();

				for (index_t i = 1; i < n; i++)
				{
					for (index_t s = 0; s < s_len; s++)
					{
						(br_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, x, noarr::lit<0>, s)) =
							(br_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, x, noarr::lit<0>, s))
							+ (diag_l | noarr::get_at<'i', 'c'>(e, i - 1, s))
								  * (br_dens_l
									 | noarr::get_at<'y', 'x', 'S', 's'>(densities, i - 1, x, noarr::lit<0>, s));
					}
				}
				for (index_t s = 0; s < s_len; s++)
				{
					(br_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, n - 1, x, noarr::lit<0>, s)) =
						(br_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, n - 1, x, noarr::lit<0>, s))
						* (diag_l | noarr::get_at<'i', 'c'>(b, n - 1, s));
				}
				for (index_t i = n - 2; i >= 0; i--)
				{
					for (index_t s = 0; s < s_len; s++)
					{
						(br_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, x, noarr::lit<0>, s)) =
							((br_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, x, noarr::lit<0>, s))
							 + c_[s]
								   * (br_dens_l
									  | noarr::get_at<'y', 'x', 'S', 's'>(densities, i + 1, x, noarr::lit<0>, s)))
							* (diag_l | noarr::get_at<'i', 'c'>(b, i, s));
					}
				}
			}
		}
	}

	// remainder
	{
		auto r_dens_l = blocked_dens_l ^ noarr::fix<'b'>(noarr::lit<1>);

		// body
		{
			auto rb_dens_l = r_dens_l ^ noarr::fix<'p'>(noarr::lit<0>);
			const index_t S_len = rb_dens_l | noarr::get_length<'S'>();
			const index_t s_len = rb_dens_l | noarr::get_length<'s'>();

#pragma omp for schedule(static) nowait
			for (index_t S = 0; S < S_len; S++)
			{
				for (index_t i = 1; i < n; i++)
				{
					for (index_t s = 0; s < s_len; s++)
					{
						(rb_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, noarr::lit<0>, S, s)) =
							(rb_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, noarr::lit<0>, S, s))
							+ (diag_l | noarr::get_at<'i', 'c'>(e, i - 1, s))
								  * (rb_dens_l
									 | noarr::get_at<'y', 'x', 'S', 's'>(densities, i - 1, noarr::lit<0>, S, s));
					}
				}
				for (index_t s = 0; s < s_len; s++)
				{
					(rb_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, n - 1, noarr::lit<0>, S, s)) =
						(rb_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, n - 1, noarr::lit<0>, S, s))
						* (diag_l | noarr::get_at<'i', 'c'>(b, n - 1, s));
				}
				for (index_t i = n - 2; i >= 0; i--)
				{
					for (index_t s = 0; s < s_len; s++)
					{
						(rb_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, noarr::lit<0>, S, s)) =
							((rb_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, noarr::lit<0>, S, s))
							 + c_[s]
								   * (rb_dens_l
									  | noarr::get_at<'y', 'x', 'S', 's'>(densities, i + 1, noarr::lit<0>, S, s)))
							* (diag_l | noarr::get_at<'i', 'c'>(b, i, s));
					}
				}
			}
		}

		// remainder
		{
			auto rr_dens_l = r_dens_l ^ noarr::fix<'p'>(noarr::lit<1>);
			const index_t s_len = rr_dens_l | noarr::get_length<'s'>();

#pragma omp single
			{
				for (index_t i = 1; i < n; i++)
				{
					for (index_t s = 0; s < s_len; s++)
					{
						(rr_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, noarr::lit<0>, noarr::lit<0>, s)) =
							(rr_dens_l
							 | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, noarr::lit<0>, noarr::lit<0>, s))
							+ (diag_l | noarr::get_at<'i', 'c'>(e, i - 1, s))
								  * (rr_dens_l
									 | noarr::get_at<'y', 'x', 'S', 's'>(densities, i - 1, noarr::lit<0>, noarr::lit<0>,
																		 s));
					}
				}
				for (index_t s = 0; s < s_len; s++)
				{
					(rr_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, n - 1, noarr::lit<0>, noarr::lit<0>, s)) =
						(rr_dens_l
						 | noarr::get_at<'y', 'x', 'S', 's'>(densities, n - 1, noarr::lit<0>, noarr::lit<0>, s))
						* (diag_l | noarr::get_at<'i', 'c'>(b, n - 1, s));
				}
				for (index_t i = n - 2; i >= 0; i--)
				{
					for (index_t s = 0; s < s_len; s++)
					{
						(rr_dens_l | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, noarr::lit<0>, noarr::lit<0>, s)) =
							((rr_dens_l
							  | noarr::get_at<'y', 'x', 'S', 's'>(densities, i, noarr::lit<0>, noarr::lit<0>, s))
							 + c_[s]
								   * (rr_dens_l
									  | noarr::get_at<'y', 'x', 'S', 's'>(densities, i + 1, noarr::lit<0>,
																		  noarr::lit<0>, s)))
							* (diag_l | noarr::get_at<'i', 'c'>(b, i, s));
					}
				}
			}
		}
	}
}

template <typename index_t, typename real_t, typename density_layout_t, typename diagonal_layout_t>
void solve_slice_y_3d(real_t* HWY_RESTRICT densities, const real_t* HWY_RESTRICT b, const real_t* HWY_RESTRICT c_,
					  const real_t* HWY_RESTRICT e, const density_layout_t dens_l, const diagonal_layout_t diag_l,
					  std::size_t s_copies, std::size_t xs_tile_size)
{
	const index_t substrate_count = dens_l | noarr::get_length<'s'>();
	const index_t n = dens_l | noarr::get_length<'y'>();
	const index_t z_len = dens_l | noarr::get_length<'z'>();

	auto blocked_dens_l = dens_l ^ noarr::merge_blocks<'x', 's', 'c'>()
						  ^ noarr::into_blocks_static<'c', 'b', 'x', 's'>(substrate_count * s_copies)
						  ^ noarr::into_blocks_static<'s', 'p', 'S', 's'>(xs_tile_size);

#pragma omp for schedule(static) nowait
	for (index_t z = 0; z < z_len; z++)
	{
		// body
		{
			auto b_dens_l = blocked_dens_l ^ noarr::fix<'b'>(noarr::lit<0>);
			const index_t x_len = b_dens_l | noarr::get_length<'x'>();

			for (index_t x = 0; x < x_len; x++)
			{
				// body
				{
					auto bb_dens_l = b_dens_l ^ noarr::fix<'p'>(noarr::lit<0>);
					const index_t S_len = bb_dens_l | noarr::get_length<'S'>();
					const index_t s_len = bb_dens_l | noarr::get_length<'s'>();

					for (index_t S = 0; S < S_len; S++)
					{
						for (index_t i = 1; i < n; i++)
						{
							for (index_t s = 0; s < s_len; s++)
							{
								(bb_dens_l | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, x, S, s)) =
									(bb_dens_l | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, x, S, s))
									+ (diag_l | noarr::get_at<'i', 'c'>(e, i - 1, s))
										  * (bb_dens_l
											 | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i - 1, x, S, s));
							}
						}
						for (index_t s = 0; s < s_len; s++)
						{
							(bb_dens_l | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, n - 1, x, S, s)) =
								(bb_dens_l | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, n - 1, x, S, s))
								* (diag_l | noarr::get_at<'i', 'c'>(b, n - 1, s));
						}
						for (index_t i = n - 2; i >= 0; i--)
						{
							for (index_t s = 0; s < s_len; s++)
							{
								(bb_dens_l | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, x, S, s)) =
									((bb_dens_l | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, x, S, s))
									 + c_[s]
										   * (bb_dens_l
											  | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i + 1, x, S, s)))
									* (diag_l | noarr::get_at<'i', 'c'>(b, i, s));
							}
						}
					}
				}

				// remainder
				{
					auto br_dens_l = b_dens_l ^ noarr::fix<'p'>(noarr::lit<1>);
					const index_t s_len = br_dens_l | noarr::get_length<'s'>();

					for (index_t i = 1; i < n; i++)
					{
						for (index_t s = 0; s < s_len; s++)
						{
							(br_dens_l | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, x, noarr::lit<0>, s)) =
								(br_dens_l
								 | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, x, noarr::lit<0>, s))
								+ (diag_l | noarr::get_at<'i', 'c'>(e, i - 1, s))
									  * (br_dens_l
										 | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i - 1, x, noarr::lit<0>,
																				  s));
						}
					}
					for (index_t s = 0; s < s_len; s++)
					{
						(br_dens_l | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, n - 1, x, noarr::lit<0>, s)) =
							(br_dens_l
							 | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, n - 1, x, noarr::lit<0>, s))
							* (diag_l | noarr::get_at<'i', 'c'>(b, n - 1, s));
					}
					for (index_t i = n - 2; i >= 0; i--)
					{
						for (index_t s = 0; s < s_len; s++)
						{
							(br_dens_l | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, x, noarr::lit<0>, s)) =
								((br_dens_l
								  | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, x, noarr::lit<0>, s))
								 + c_[s]
									   * (br_dens_l
										  | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i + 1, x,
																				   noarr::lit<0>, s)))
								* (diag_l | noarr::get_at<'i', 'c'>(b, i, s));
						}
					}
				}
			}
		}

		// remainder
		{
			auto r_dens_l = blocked_dens_l ^ noarr::fix<'b'>(noarr::lit<1>);

			// body
			{
				auto rb_dens_l = r_dens_l ^ noarr::fix<'p'>(noarr::lit<0>);
				const index_t S_len = rb_dens_l | noarr::get_length<'S'>();
				const index_t s_len = rb_dens_l | noarr::get_length<'s'>();

				for (index_t S = 0; S < S_len; S++)
				{
					for (index_t i = 1; i < n; i++)
					{
						for (index_t s = 0; s < s_len; s++)
						{
							(rb_dens_l | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, noarr::lit<0>, S, s)) =
								(rb_dens_l
								 | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, noarr::lit<0>, S, s))
								+ (diag_l | noarr::get_at<'i', 'c'>(e, i - 1, s))
									  * (rb_dens_l
										 | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i - 1, noarr::lit<0>, S,
																				  s));
						}
					}
					for (index_t s = 0; s < s_len; s++)
					{
						(rb_dens_l | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, n - 1, noarr::lit<0>, S, s)) =
							(rb_dens_l
							 | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, n - 1, noarr::lit<0>, S, s))
							* (diag_l | noarr::get_at<'i', 'c'>(b, n - 1, s));
					}
					for (index_t i = n - 2; i >= 0; i--)
					{
						for (index_t s = 0; s < s_len; s++)
						{
							(rb_dens_l | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, noarr::lit<0>, S, s)) =
								((rb_dens_l
								  | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, noarr::lit<0>, S, s))
								 + c_[s]
									   * (rb_dens_l
										  | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i + 1, noarr::lit<0>,
																				   S, s)))
								* (diag_l | noarr::get_at<'i', 'c'>(b, i, s));
						}
					}
				}
			}

			// remainder
			{
				auto rr_dens_l = r_dens_l ^ noarr::fix<'p'>(noarr::lit<1>);
				const index_t s_len = rr_dens_l | noarr::get_length<'s'>();

				for (index_t i = 1; i < n; i++)
				{
					for (index_t s = 0; s < s_len; s++)
					{
						(rr_dens_l
						 | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, noarr::lit<0>, noarr::lit<0>, s)) =
							(rr_dens_l
							 | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, noarr::lit<0>, noarr::lit<0>, s))
							+ (diag_l | noarr::get_at<'i', 'c'>(e, i - 1, s))
								  * (rr_dens_l
									 | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i - 1, noarr::lit<0>,
																			  noarr::lit<0>, s));
					}
				}
				for (index_t s = 0; s < s_len; s++)
				{
					(rr_dens_l
					 | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, n - 1, noarr::lit<0>, noarr::lit<0>, s)) =
						(rr_dens_l
						 | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, n - 1, noarr::lit<0>, noarr::lit<0>, s))
						* (diag_l | noarr::get_at<'i', 'c'>(b, n - 1, s));
				}
				for (index_t i = n - 2; i >= 0; i--)
				{
					for (index_t s = 0; s < s_len; s++)
					{
						(rr_dens_l
						 | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, noarr::lit<0>, noarr::lit<0>, s)) =
							((rr_dens_l
							  | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i, noarr::lit<0>, noarr::lit<0>,
																	   s))
							 + c_[s]
								   * (rr_dens_l
									  | noarr::get_at<'z', 'y', 'x', 'S', 's'>(densities, z, i + 1, noarr::lit<0>,
																			   noarr::lit<0>, s)))
							* (diag_l | noarr::get_at<'i', 'c'>(b, i, s));
					}
				}
			}
		}
	}
}

template <typename index_t, typename real_t, typename density_layout_t, typename diagonal_layout_t>
void solve_slice_z_3d(real_t* HWY_RESTRICT densities, const real_t* HWY_RESTRICT b, const real_t* HWY_RESTRICT c_,
					  const real_t* HWY_RESTRICT e, const density_layout_t dens_l, const diagonal_layout_t diag_l,
					  std::size_t s_copies, std::size_t xs_tile_size)
{
	const index_t substrate_count = dens_l | noarr::get_length<'s'>();
	const index_t n = dens_l | noarr::get_length<'z'>();
	const index_t y_len = dens_l | noarr::get_length<'y'>();

	auto blocked_dens_l = dens_l ^ noarr::merge_blocks<'x', 's', 'c'>()
						  ^ noarr::into_blocks_static<'c', 'b', 'x', 's'>(substrate_count * s_copies)
						  ^ noarr::into_blocks_static<'s', 'p', 'S', 's'>(xs_tile_size);

#pragma omp for schedule(static) nowait
	for (index_t y = 0; y < y_len; y++)
	{
		// body
		{
			auto b_dens_l = blocked_dens_l ^ noarr::fix<'b'>(noarr::lit<0>);
			const index_t x_len = b_dens_l | noarr::get_length<'x'>();

			for (index_t x = 0; x < x_len; x++)
			{
				// body
				{
					auto bb_dens_l = b_dens_l ^ noarr::fix<'p'>(noarr::lit<0>);
					const index_t S_len = bb_dens_l | noarr::get_length<'S'>();
					const index_t s_len = bb_dens_l | noarr::get_length<'s'>();

					for (index_t S = 0; S < S_len; S++)
					{
						for (index_t i = 1; i < n; i++)
						{
							for (index_t s = 0; s < s_len; s++)
							{
								(bb_dens_l | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, x, S, s)) =
									(bb_dens_l | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, x, S, s))
									+ (diag_l | noarr::get_at<'i', 'c'>(e, i - 1, s))
										  * (bb_dens_l
											 | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i - 1, x, S, s));
							}
						}
						for (index_t s = 0; s < s_len; s++)
						{
							(bb_dens_l | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, n - 1, x, S, s)) =
								(bb_dens_l | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, n - 1, x, S, s))
								* (diag_l | noarr::get_at<'i', 'c'>(b, n - 1, s));
						}
						for (index_t i = n - 2; i >= 0; i--)
						{
							for (index_t s = 0; s < s_len; s++)
							{
								(bb_dens_l | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, x, S, s)) =
									((bb_dens_l | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, x, S, s))
									 + c_[s]
										   * (bb_dens_l
											  | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i + 1, x, S, s)))
									* (diag_l | noarr::get_at<'i', 'c'>(b, i, s));
							}
						}
					}
				}

				// remainder
				{
					auto br_dens_l = b_dens_l ^ noarr::fix<'p'>(noarr::lit<1>);
					const index_t s_len = br_dens_l | noarr::get_length<'s'>();

					for (index_t i = 1; i < n; i++)
					{
						for (index_t s = 0; s < s_len; s++)
						{
							(br_dens_l | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, x, noarr::lit<0>, s)) =
								(br_dens_l
								 | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, x, noarr::lit<0>, s))
								+ (diag_l | noarr::get_at<'i', 'c'>(e, i - 1, s))
									  * (br_dens_l
										 | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i - 1, x, noarr::lit<0>,
																				  s));
						}
					}
					for (index_t s = 0; s < s_len; s++)
					{
						(br_dens_l | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, n - 1, x, noarr::lit<0>, s)) =
							(br_dens_l
							 | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, n - 1, x, noarr::lit<0>, s))
							* (diag_l | noarr::get_at<'i', 'c'>(b, n - 1, s));
					}
					for (index_t i = n - 2; i >= 0; i--)
					{
						for (index_t s = 0; s < s_len; s++)
						{
							(br_dens_l | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, x, noarr::lit<0>, s)) =
								((br_dens_l
								  | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, x, noarr::lit<0>, s))
								 + c_[s]
									   * (br_dens_l
										  | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i + 1, x,
																				   noarr::lit<0>, s)))
								* (diag_l | noarr::get_at<'i', 'c'>(b, i, s));
						}
					}
				}
			}
		}

		// remainder
		{
			auto r_dens_l = blocked_dens_l ^ noarr::fix<'b'>(noarr::lit<1>);

			// body
			{
				auto rb_dens_l = r_dens_l ^ noarr::fix<'p'>(noarr::lit<0>);
				const index_t S_len = rb_dens_l | noarr::get_length<'S'>();
				const index_t s_len = rb_dens_l | noarr::get_length<'s'>();

				for (index_t S = 0; S < S_len; S++)
				{
					for (index_t i = 1; i < n; i++)
					{
						for (index_t s = 0; s < s_len; s++)
						{
							(rb_dens_l | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, noarr::lit<0>, S, s)) =
								(rb_dens_l
								 | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, noarr::lit<0>, S, s))
								+ (diag_l | noarr::get_at<'i', 'c'>(e, i - 1, s))
									  * (rb_dens_l
										 | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i - 1, noarr::lit<0>, S,
																				  s));
						}
					}
					for (index_t s = 0; s < s_len; s++)
					{
						(rb_dens_l | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, n - 1, noarr::lit<0>, S, s)) =
							(rb_dens_l
							 | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, n - 1, noarr::lit<0>, S, s))
							* (diag_l | noarr::get_at<'i', 'c'>(b, n - 1, s));
					}
					for (index_t i = n - 2; i >= 0; i--)
					{
						for (index_t s = 0; s < s_len; s++)
						{
							(rb_dens_l | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, noarr::lit<0>, S, s)) =
								((rb_dens_l
								  | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, noarr::lit<0>, S, s))
								 + c_[s]
									   * (rb_dens_l
										  | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i + 1, noarr::lit<0>,
																				   S, s)))
								* (diag_l | noarr::get_at<'i', 'c'>(b, i, s));
						}
					}
				}
			}

			// remainder
			{
				auto rr_dens_l = r_dens_l ^ noarr::fix<'p'>(noarr::lit<1>);
				const index_t s_len = rr_dens_l | noarr::get_length<'s'>();

				for (index_t i = 1; i < n; i++)
				{
					for (index_t s = 0; s < s_len; s++)
					{
						(rr_dens_l
						 | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, noarr::lit<0>, noarr::lit<0>, s)) =
							(rr_dens_l
							 | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, noarr::lit<0>, noarr::lit<0>, s))
							+ (diag_l | noarr::get_at<'i', 'c'>(e, i - 1, s))
								  * (rr_dens_l
									 | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i - 1, noarr::lit<0>,
																			  noarr::lit<0>, s));
					}
				}
				for (index_t s = 0; s < s_len; s++)
				{
					(rr_dens_l
					 | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, n - 1, noarr::lit<0>, noarr::lit<0>, s)) =
						(rr_dens_l
						 | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, n - 1, noarr::lit<0>, noarr::lit<0>, s))
						* (diag_l | noarr::get_at<'i', 'c'>(b, n - 1, s));
				}
				for (index_t i = n - 2; i >= 0; i--)
				{
					for (index_t s = 0; s < s_len; s++)
					{
						(rr_dens_l
						 | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, noarr::lit<0>, noarr::lit<0>, s)) =
							((rr_dens_l
							  | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i, noarr::lit<0>, noarr::lit<0>,
																	   s))
							 + c_[s]
								   * (rr_dens_l
									  | noarr::get_at<'y', 'z', 'x', 'S', 's'>(densities, y, i + 1, noarr::lit<0>,
																			   noarr::lit<0>, s)))
							* (diag_l | noarr::get_at<'i', 'c'>(b, i, s));
					}
				}
			}
		}
	}
}
} // namespace

void diffusion_solver::solve()
{
	if (problem.dims == 1)
	{
		for (index_t i = 0; i < problem.iterations; i++)
		{
			solve_slice_x_1d<sindex_t>(this->substrates_.get(), bx_.get(), cx_.get(), ex_.get(),
									   get_substrates_layout<1>(), get_diagonal_layout(problem, problem.nx));
#pragma omp barrier
		}
	}
	else if (problem.dims == 2)
	{
		for (index_t i = 0; i < problem.iterations; i++)
		{
			solve_slice_x_2d_and_3d<sindex_t>(this->substrates_.get(), bx_.get(), cx_.get(), ex_.get(),
											  get_substrates_layout<2>() ^ noarr::rename<'y', 'm'>(),
											  get_diagonal_layout(problem, problem.nx));
#pragma omp barrier
			solve_slice_y_2d<sindex_t>(this->substrates_.get(), by_.get(), cy_.get(), ey_.get(),
									   get_substrates_layout<2>(),
									   get_diagonal_layout_c(problem, problem.ny, (index_t)substrate_copies_),
									   substrate_copies_, xs_tile_size_);
#pragma omp barrier
		}
	}
	else if (problem.dims == 3)
	{
		for (index_t i = 0; i < problem.iterations; i++)
		{
			solve_slice_x_2d_and_3d<sindex_t>(this->substrates_.get(), bx_.get(), cx_.get(), ex_.get(),
											  get_substrates_layout<3>() ^ noarr::merge_blocks<'z', 'y', 'm'>(),
											  get_diagonal_layout(problem, problem.nx));
#pragma omp barrier
			solve_slice_y_3d<sindex_t>(this->substrates_.get(), by_.get(), cy_.get(), ey_.get(),
									   get_substrates_layout<3>(),
									   get_diagonal_layout_c(problem, problem.ny, (index_t)substrate_copies_),
									   substrate_copies_, xs_tile_size_);
#pragma omp barrier
			solve_slice_z_3d<sindex_t>(this->substrates_.get(), bz_.get(), cz_.get(), ez_.get(),
									   get_substrates_layout<3>(),
									   get_diagonal_layout_c(problem, problem.ny, (index_t)substrate_copies_),
									   substrate_copies_, xs_tile_size_);
#pragma omp barrier
		}
	}
}
