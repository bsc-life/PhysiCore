#pragma once

#include <memory>

#include <common/types.h>
#include <hwy/aligned_allocator.h>
#include <noarr/structures_extended.hpp>

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

namespace physicore::biofvm::kernels::openmp_solver {

class position_solver
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
	void prepare(const microenvironment& m, index_t iterations);

	void initialize();

	void solve();

	static void update_cell_forces(environment& e);

	static void update_cell_neighbors(environment& e);

	static void update_motility(environment& e);

	static void update_basement_membrane_interactions(environment& e);

	static void update_spring_attachments(environment& e);

	static void update_positions(environment& e);
};

} // namespace physicore::biofvm::kernels::openmp_solver
