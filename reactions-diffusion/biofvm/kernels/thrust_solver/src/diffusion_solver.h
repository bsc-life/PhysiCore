#pragma once

#include <biofvm/microenvironment.h>
#include <thrust/device_ptr.h>

#include "data_manager.h"

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
*/

namespace physicore::biofvm::kernels::thrust_solver {

class diffusion_solver
{
	thrust::device_ptr<real_t> bx_, cx_, ex_;
	thrust::device_ptr<real_t> by_, cy_, ey_;
	thrust::device_ptr<real_t> bz_, cz_, ez_;

	index_t substrate_factor_ = 0;

	index_t ns_ = 0;
	index_t nx_ = 0;
	index_t ny_ = 0;
	index_t nz_ = 0;

	thrust::device_ptr<real_t> substrate_densities_;

	thrust::device_ptr<real_t> initial_conditions_;

public:
	diffusion_solver() = default;
	diffusion_solver(const diffusion_solver&) = delete;
	diffusion_solver& operator=(const diffusion_solver&) = delete;

	static void precompute_values(thrust::device_ptr<real_t>& b, thrust::device_ptr<real_t>& c,
								  thrust::device_ptr<real_t>& e, index_t shape, index_t dims, index_t n,
								  const microenvironment& m, index_t copies);

	template <std::size_t dims = 3>
	auto get_substrates_layout() const
	{
		if constexpr (dims == 1)
			return noarr::scalar<real_t>() ^ noarr::vectors<'s', 'x'>(ns_, nx_);
		else if constexpr (dims == 2)
			return noarr::scalar<real_t>() ^ noarr::vectors<'s', 'x', 'y'>(ns_, nx_, ny_);
		else if constexpr (dims == 3)
			return noarr::scalar<real_t>() ^ noarr::vectors<'s', 'x', 'y', 'z'>(ns_, nx_, ny_, nz_);
	}

	thrust::device_ptr<real_t> get_substrates_pointer();

	void initialize(microenvironment& m);
	void initialize(microenvironment& m, index_t substrate_factor);

	void deinitialize();

	void solve();

	~diffusion_solver();
};

} // namespace physicore::biofvm::kernels::thrust_solver
