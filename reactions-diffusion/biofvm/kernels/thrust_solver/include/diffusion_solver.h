#pragma once

#include "../../../include/microenvironment.h"
#include "device_manager.h"

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
b_1'  == 1/b_1
b_i'  == 1/(b_i - a_i*c_i*b_(i-1)')                           1 <  i <= n
e_i   == a_i*b_(i-1)'                                         1 <  i <= n

Then, the forward substitution is as follows (n multiplication + n subtractions):
d_i'  == d_i - e_i*d_(i-1)                                    1 <  i <= n
The backpropagation (2n multiplication + n subtractions):
d_n'' == d_n'/b_n'
d_i'' == (d_i' - c_i*d_(i+1)'')*b_i'                          n >  i >= 1

Optimizations:
- The memory accesses of all dimension slices are continuous for the least amount of cache misses. This is done by
reordering the loops for the thomas solver of y and z axes.
- Better vectorization for small enough number of substrates. Since compiler generated vectorization is done over the
substrates, when we have 1 or 2 substrates, vector capability is not used to the full. We can make it better for y and z
axes by merging loops over substrates (s) and x dimension. When we merge these loops, the vectorization can be done over
the whole 'plane' s*x rather than just s. For it to work properly, helper vectors b, c, e had to be modified such that
their data span over x*s, not only over s (so they got copied substrate_factor_ times in s dimension).
*/

namespace physicore::biofvm::kernels::thrust_solver {

class diffusion_solver
{
	device_vector<real_t> bx_, cx_, ex_;
	device_vector<real_t> by_, cy_, ey_;
	device_vector<real_t> bz_, cz_, ez_;

	index_t substrate_factor_;

	index_t ns_;
	index_t nx_;
	index_t ny_;
	index_t nz_;

	device_vector<real_t> substrate_densities_;

public:
	static void precompute_values(device_vector<real_t>& b, device_vector<real_t>& c, device_vector<real_t>& e,
								  index_t shape, index_t dims, index_t n, const microenvironment& m, index_t copies);

	template <std::size_t dims = 3>
	auto get_substrates_layout()
	{
		if constexpr (dims == 1)
			return noarr::scalar<real_t>() ^ noarr::vectors<'s', 'x'>(ns_, nx_);
		else if constexpr (dims == 2)
			return noarr::scalar<real_t>() ^ noarr::vectors<'s', 'x', 'y'>(ns_, nx_, ny_);
		else if constexpr (dims == 3)
			return noarr::scalar<real_t>() ^ noarr::vectors<'s', 'x', 'y', 'z'>(ns_, nx_, ny_, nz_);
	}

	real_t* get_substrates_pointer();

	void initialize(microenvironment& m);
	void initialize(microenvironment& m, index_t substrate_factor);

	void solve();
};

} // namespace physicore::biofvm::kernels::thrust_solver
