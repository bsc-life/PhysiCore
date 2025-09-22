#pragma once

#include <memory>

#include <hwy/aligned_allocator.h>
#include <noarr/structures_extended.hpp>

#include "types.h"
#include "problem.h"

/*
The diffusion is the problem of solving tridiagonal matrix system with these coeficients:
For dimension x:
a_i  == -dt*diffusion_coefs/dx^2                              1 <= i <= n
b_1  == 1 + dt*decay_rates/dims + dt*diffusion_coefs/dx^2
b_i  == 1 + dt*decay_rates/dims + 2*dt*diffusion_coefs/dx^2   1 <  i <  n
b_n  == 1 + dt*decay_rates/dims + dt*diffusion_coefs/dx^2
c_i  == -dt*diffusion_coefs/dx^2                              1 <= i <= n
d_i  == current diffusion rates
For dimension y/z (if they exist):
substitute dx accordingly to dy/dz

Since the matrix is constant for multiple right hand sides, we precompute its values in the following way:
a_i' = c_i' = -a_i
b_1'  == 1/b_1
b_i'  == 1/(b_i - a_i'*c_i'*b_(i-1)')                         1 <  i <= n
e_i   == a_i'*b_(i-1)'                                        1 <  i <= n

Then, the forward substitution is as follows (n FMAs):
d_i'  == d_i + e_i*d_(i-1)                                    1 <  i <= n
The backpropagation (n multiplications + n FMAs):
d_n'' == d_n'/b_n'
d_i'' == (d_i' + c_i*d_(i+1)'')*b_i'                          n >  i >= 1

Optimizations:
- Precomputed a_i, b_i', e_i
- Aligned memory for x dimension (tunable by 'alignment_size')
- Better temporal locality of memory accesses - sx plane is divided into smaller tiles (tunable by 'xs_tile_size') and
y/z dimensions are solved alongside tiled xs dimension
*/

namespace physicore::biofvm::kernels::cpu {

class diffusion_solver
{
private:
	problem_t problem;

	std::unique_ptr<real_t[]> bx_, cx_, ex_;
	std::unique_ptr<real_t[]> by_, cy_, ey_;
	std::unique_ptr<real_t[]> bz_, cz_, ez_;

	std::size_t xs_tile_size_ = 48;
	std::size_t alignment_size_ = HWY_ALIGNMENT;

	std::size_t substrate_copies_;

	hwy::AlignedUniquePtr<real_t[]> substrates_;

	void precompute_values(std::unique_ptr<real_t[]>& b, std::unique_ptr<real_t[]>& c, std::unique_ptr<real_t[]>& e,
						   index_t shape, index_t dims, index_t n, index_t copies);

public:
	template <std::size_t dims = 3>
	auto get_substrates_layout()
	{
		std::size_t xs_size = problem.nx * problem.substrates_count * sizeof(real_t);
		std::size_t xs_size_padded = (xs_size + alignment_size_ - 1) / alignment_size_ * alignment_size_;
		xs_size_padded /= sizeof(real_t);

		if constexpr (dims == 1)
			return noarr::scalar<real_t>() ^ noarr::vectors<'x'>(xs_size_padded)
				   ^ noarr::into_blocks_static<'x', 'b', 'x', 's'>(problem.substrates_count)
				   ^ noarr::fix<'b'>(noarr::lit<0>) ^ noarr::slice<'x'>(problem.nx);
		else if constexpr (dims == 2)
			return noarr::scalar<real_t>() ^ noarr::vectors<'x', 'y'>(xs_size_padded, problem.ny)
				   ^ noarr::into_blocks_static<'x', 'b', 'x', 's'>(problem.substrates_count)
				   ^ noarr::fix<'b'>(noarr::lit<0>) ^ noarr::slice<'x'>(problem.nx);
		else if constexpr (dims == 3)
			return noarr::scalar<real_t>() ^ noarr::vectors<'x', 'y', 'z'>(xs_size_padded, problem.ny, problem.nz)
				   ^ noarr::into_blocks_static<'x', 'b', 'x', 's'>(problem.substrates_count)
				   ^ noarr::fix<'b'>(noarr::lit<0>) ^ noarr::slice<'x'>(problem.nx);
	}

	real_t* get_substrates_pointer();

	void prepare(const microenvironment& m, index_t iterations);

	void initialize();

	void solve();
};

} // namespace physicore::biofvm::kernels::cpu
