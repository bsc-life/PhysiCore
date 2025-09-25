#include "cell_solver.h"

#include <limits>

#include <noarr/structures_extended.hpp>

using namespace physicore;
using namespace physicore::biofvm;
using namespace physicore::biofvm::kernels::cpu;


static constexpr index_t no_ballot = std::numeric_limits<index_t>::max();

template <index_t dims>
auto fix_dims(const real_t* cell_position, const cartesian_mesh& m)
{
	std::array<index_t, 3> voxel_index = m.voxel_position(std::span<const real_t, dims>(cell_position, dims));
	return noarr::fix<'x'>(voxel_index[0]) ^ noarr::fix<'y'>(voxel_index[1]) ^ noarr::fix<'z'>(voxel_index[2]);
}

template <index_t dims>
void clear_ballots(const auto ballot_l, const real_t* HWY_RESTRICT cell_positions,
				   std::atomic<index_t>* HWY_RESTRICT ballots, std::atomic<real_t>* HWY_RESTRICT reduced_numerators,
				   std::atomic<real_t>* HWY_RESTRICT reduced_denominators,
				   std::atomic<real_t>* HWY_RESTRICT reduced_factors, index_t n, const cartesian_mesh& m,
				   index_t substrate_densities)
{
#pragma omp for
	for (index_t i = 0; i < n; i++)
	{
		auto b_l = ballot_l ^ fix_dims<dims>(cell_positions + dims * i, m);

		auto& b = b_l | noarr::get_at(ballots);

		b.store(no_ballot, std::memory_order_relaxed);

		for (index_t s = 0; s < substrate_densities; s++)
		{
			reduced_numerators[i * substrate_densities + s].store(0, std::memory_order_relaxed);
			reduced_denominators[i * substrate_densities + s].store(0, std::memory_order_relaxed);
			reduced_factors[i * substrate_densities + s].store(0, std::memory_order_relaxed);
		}
	}
}

void compute_intermediates(real_t* HWY_RESTRICT numerators, real_t* HWY_RESTRICT denominators,
						   real_t* HWY_RESTRICT factors, const real_t* HWY_RESTRICT secretion_rates,
						   const real_t* HWY_RESTRICT uptake_rates, const real_t* HWY_RESTRICT saturation_densities,
						   const real_t* HWY_RESTRICT net_export_rates, const real_t* HWY_RESTRICT cell_volumes,
						   real_t voxel_volume, real_t time_step, index_t n, index_t substrates_count)
{
#pragma omp for
	for (index_t i = 0; i < n; i++)
	{
		for (index_t s = 0; s < substrates_count; s++)
		{
			numerators[i * substrates_count + s] = secretion_rates[i * substrates_count + s]
												   * saturation_densities[i * substrates_count + s] * time_step
												   * cell_volumes[i] / voxel_volume;

			denominators[i * substrates_count + s] =
				(uptake_rates[i * substrates_count + s] + secretion_rates[i * substrates_count + s]) * time_step
				* cell_volumes[i] / voxel_volume;

			factors[i * substrates_count + s] = net_export_rates[i * substrates_count + s] * time_step / voxel_volume;
		}
	}
}

template <index_t dims>
void ballot_and_sum(const auto ballot_l, std::atomic<real_t>* HWY_RESTRICT reduced_numerators,
					std::atomic<real_t>* HWY_RESTRICT reduced_denominators,
					std::atomic<real_t>* HWY_RESTRICT reduced_factors, const real_t* HWY_RESTRICT numerators,
					const real_t* HWY_RESTRICT denominators, const real_t* HWY_RESTRICT factors,
					const real_t* HWY_RESTRICT cell_positions, std::atomic<index_t>* HWY_RESTRICT ballots, index_t n,
					index_t substrates_count, const cartesian_mesh& m, std::atomic<bool>* HWY_RESTRICT is_conflict)
{
#pragma omp for
	for (index_t i = 0; i < n; i++)
	{
		auto b_l = ballot_l ^ fix_dims<dims>(cell_positions + dims * i, m);

		auto& b = b_l | noarr::get_at(ballots);

		auto expected = no_ballot;
		bool success = b.compare_exchange_strong(expected, i, std::memory_order_acq_rel, std::memory_order_acquire);

		if (success)
		{
			for (index_t s = 0; s < substrates_count; s++)
			{
				reduced_numerators[i * substrates_count + s].fetch_add(numerators[i * substrates_count + s],
																	   std::memory_order_relaxed);
				reduced_denominators[i * substrates_count + s].fetch_add(denominators[i * substrates_count + s] + 1,
																		 std::memory_order_relaxed);
				reduced_factors[i * substrates_count + s].fetch_add(factors[i * substrates_count + s],
																	std::memory_order_relaxed);
			}
		}
		else
		{
			is_conflict[0].store(true, std::memory_order_relaxed);

			for (index_t s = 0; s < substrates_count; s++)
			{
				reduced_numerators[expected * substrates_count + s].fetch_add(numerators[i * substrates_count + s],
																			  std::memory_order_relaxed);
				reduced_denominators[expected * substrates_count + s].fetch_add(denominators[i * substrates_count + s],
																				std::memory_order_relaxed);
				reduced_factors[expected * substrates_count + s].fetch_add(factors[i * substrates_count + s],
																		   std::memory_order_relaxed);
			}
		}
	}
}

template <typename density_layout_t>
void compute_internalized(real_t* HWY_RESTRICT internalized_substrates, const real_t* HWY_RESTRICT substrate_densities,
						  const real_t* HWY_RESTRICT numerator, const real_t* HWY_RESTRICT denominator,
						  const real_t* HWY_RESTRICT factor, real_t voxel_volume, density_layout_t dens_l)
{
	const index_t substrates_count = dens_l | noarr::get_length<'s'>();

	for (index_t s = 0; s < substrates_count; s++)
	{
		internalized_substrates[s] -=
			voxel_volume
			* (((-denominator[s]) * (dens_l | noarr::get_at<'s'>(substrate_densities, s)) + numerator[s])
				   / (1 + denominator[s])
			   + factor[s]);
	}
}

template <typename density_layout_t>
void compute_densities(real_t* HWY_RESTRICT substrate_densities, const std::atomic<real_t>* HWY_RESTRICT numerator,
					   const std::atomic<real_t>* HWY_RESTRICT denominator,
					   const std::atomic<real_t>* HWY_RESTRICT factor, bool has_ballot, density_layout_t dens_l)
{
	const index_t substrates_count = dens_l | noarr::get_length<'s'>();

	if (has_ballot)
	{
		for (index_t s = 0; s < substrates_count; s++)
		{
			(dens_l | noarr::get_at<'s'>(substrate_densities, s)) =
				((dens_l | noarr::get_at<'s'>(substrate_densities, s)) + numerator[s].load(std::memory_order_relaxed))
					/ denominator[s].load(std::memory_order_relaxed)
				+ factor[s].load(std::memory_order_relaxed);
		}
	}
}

template <typename density_layout_t>
void compute_fused(real_t* HWY_RESTRICT substrate_densities, real_t* HWY_RESTRICT internalized_substrates,
				   const std::atomic<real_t>* HWY_RESTRICT numerator,
				   const std::atomic<real_t>* HWY_RESTRICT denominator, const std::atomic<real_t>* HWY_RESTRICT factor,
				   real_t voxel_volume, density_layout_t dens_l)
{
	const index_t substrates_count = dens_l | noarr::get_length<'s'>();

	for (index_t s = 0; s < substrates_count; s++)
	{
		internalized_substrates[s] -= voxel_volume
									  * (((1 - denominator[s].load(std::memory_order_relaxed))
											  * (dens_l | noarr::get_at<'s'>(substrate_densities, s))
										  + numerator[s].load(std::memory_order_relaxed))
											 / (denominator[s].load(std::memory_order_relaxed))
										 + factor[s].load(std::memory_order_relaxed));

		(dens_l | noarr::get_at<'s'>(substrate_densities, s)) =
			((dens_l | noarr::get_at<'s'>(substrate_densities, s)) + numerator[s].load(std::memory_order_relaxed))
				/ denominator[s].load(std::memory_order_relaxed)
			+ factor[s].load(std::memory_order_relaxed);
	}
}

template <index_t dims>
void compute_result(const auto dens_l, const auto ballot_l, agent_data& data, const cartesian_mesh& mesh,
					real_t* substrates, const std::atomic<real_t>* reduced_numerators,
					const std::atomic<real_t>* reduced_denominators, const std::atomic<real_t>* reduced_factors,
					const real_t* numerators, const real_t* denominators, const real_t* factors,
					const std::atomic<index_t>* ballots, bool with_internalized, bool is_conflict)
{
	auto voxel_volume = (real_t)mesh.voxel_volume(); // expecting that voxel volume is the same for all voxels

	if (with_internalized && !is_conflict)
	{
#pragma omp for
		for (index_t i = 0; i < data.base_data.agents_count; i++)
		{
			auto fixed_dims = fix_dims<dims>(data.base_data.positions.data() + i * dims, mesh);

			compute_fused(substrates, data.internalized_substrates.data() + i * data.substrate_count,
						  reduced_numerators + i * data.substrate_count,
						  reduced_denominators + i * data.substrate_count, reduced_factors + i * data.substrate_count,
						  voxel_volume, dens_l ^ fixed_dims);
		}

		return;
	}

	if (with_internalized)
	{
#pragma omp for
		for (index_t i = 0; i < data.base_data.agents_count; i++)
		{
			auto fixed_dims = fix_dims<dims>(data.base_data.positions.data() + i * dims, mesh);

			compute_internalized(data.internalized_substrates.data() + i * data.substrate_count, substrates,
								 numerators + i * data.substrate_count, denominators + i * data.substrate_count,
								 factors + i * data.substrate_count, voxel_volume, dens_l ^ fixed_dims);
		}
	}

#pragma omp for
	for (index_t i = 0; i < data.base_data.agents_count; i++)
	{
		auto fixed_dims = fix_dims<dims>(data.base_data.positions.data() + i * dims, mesh);

		auto ballot = ((ballot_l ^ fixed_dims) | noarr::get_at(ballots)).load(std::memory_order_relaxed);
		compute_densities(substrates, reduced_numerators + i * data.substrate_count,
						  reduced_denominators + i * data.substrate_count, reduced_factors + i * data.substrate_count,
						  ballot == i, dens_l ^ fixed_dims);
	}
}

template <index_t dims>
void simulate(const auto dens_l, const auto ballot_l, agent_data& data, microenvironment& m, real_t* substrates,
			  std::atomic<real_t>* reduced_numerators, std::atomic<real_t>* reduced_denominators,
			  std::atomic<real_t>* reduced_factors, real_t* numerators, real_t* denominators, real_t* factors,
			  std::atomic<index_t>* ballots, bool recompute, bool with_internalized,
			  std::atomic<bool>* HWY_RESTRICT is_conflict)
{
	if (recompute)
	{
		compute_intermediates(numerators, denominators, factors, data.secretion_rates.data(), data.uptake_rates.data(),
							  data.saturation_densities.data(), data.net_export_rates.data(), data.volumes.data(),
							  (real_t)m.mesh.voxel_volume(), m.diffusion_timestep, data.base_data.agents_count,
							  data.substrate_count);

		clear_ballots<dims>(ballot_l, data.base_data.positions.data(), ballots, reduced_numerators,
							reduced_denominators, reduced_factors, data.base_data.agents_count, m.mesh,
							data.substrate_count);

		ballot_and_sum<dims>(ballot_l, reduced_numerators, reduced_denominators, reduced_factors, numerators,
							 denominators, factors, data.base_data.positions.data(), ballots,
							 data.base_data.agents_count, data.substrate_count, m.mesh, is_conflict);
	}

	compute_result<dims>(dens_l, ballot_l, data, m.mesh, substrates, reduced_numerators, reduced_denominators,
						 reduced_factors, numerators, denominators, factors, ballots, with_internalized,
						 is_conflict[0].load(std::memory_order_relaxed));
}

void cell_solver::simulate_secretion_and_uptake(microenvironment& m, diffusion_solver& d_solver, bool recompute)
{
	real_t* substrates = d_solver.get_substrates_pointer();

#pragma omp single
	if (recompute)
	{
		resize(m);
		is_conflict_.store(false, std::memory_order_relaxed);
	}

	switch (m.mesh.dims)
	{
		case 1: {
			const auto dens_l = d_solver.get_substrates_layout<1>();
			const auto ballot_l = noarr::scalar<std::atomic<index_t>>() ^ noarr::vectors<'x'>(m.mesh.grid_shape[0]);

			simulate<1>(dens_l, ballot_l, retrieve_agent_data(m.agents), m, substrates, reduced_numerators_.get(),
						reduced_denominators_.get(), reduced_factors_.get(), numerators_.data(), denominators_.data(),
						factors_.data(), ballots_.get(), recompute, compute_internalized_substrates_, &is_conflict_);
			return;
		}
		case 2: {
			const auto dens_l = d_solver.get_substrates_layout<2>();
			const auto ballot_l = noarr::scalar<std::atomic<index_t>>()
								  ^ noarr::vectors<'x', 'y'>(m.mesh.grid_shape[0], m.mesh.grid_shape[1]);

			simulate<2>(dens_l, ballot_l, retrieve_agent_data(m.agents), m, substrates, reduced_numerators_.get(),
						reduced_denominators_.get(), reduced_factors_.get(), numerators_.data(), denominators_.data(),
						factors_.data(), ballots_.get(), recompute, compute_internalized_substrates_, &is_conflict_);
			return;
		}
		case 3: {
			const auto dens_l = d_solver.get_substrates_layout<3>();
			const auto ballot_l =
				noarr::scalar<std::atomic<index_t>>()
				^ noarr::vectors<'x', 'y', 'z'>(m.mesh.grid_shape[0], m.mesh.grid_shape[1], m.mesh.grid_shape[2]);

			simulate<3>(dens_l, ballot_l, retrieve_agent_data(m.agents), m, substrates, reduced_numerators_.get(),
						reduced_denominators_.get(), reduced_factors_.get(), numerators_.data(), denominators_.data(),
						factors_.data(), ballots_.get(), recompute, compute_internalized_substrates_, &is_conflict_);
			return;
		}
		default:
			assert(false);
			return;
	}
}

template <typename density_layout_t>
void release_internal(real_t* HWY_RESTRICT substrate_densities, real_t* HWY_RESTRICT internalized_substrates,
					  const real_t* HWY_RESTRICT fraction_released_at_death, real_t voxel_volume,
					  density_layout_t dens_l)
{
	const index_t substrates_count = dens_l | noarr::get_length<'s'>();

	for (index_t s = 0; s < substrates_count; s++)
	{
		std::atomic_ref<real_t>(dens_l | noarr::get_at<'s'>(substrate_densities, s))
			.fetch_add(internalized_substrates[s] * fraction_released_at_death[s] / voxel_volume,
					   std::memory_order_relaxed);

		internalized_substrates[s] = 0;
	}
}

template <index_t dims>
void release_dim(const auto dens_l, agent_data& data, const cartesian_mesh& mesh, real_t* substrates, index_t index)
{
	auto voxel_volume = (real_t)mesh.voxel_volume(); // expecting that voxel volume is the same for all voxels

	release_internal(substrates, data.internalized_substrates.data() + index * data.substrate_count,
					 data.fraction_released_at_death.data() + index * data.substrate_count, voxel_volume,
					 dens_l ^ fix_dims<dims>(data.base_data.positions.data() + index * dims, mesh));
}

void cell_solver::release_internalized_substrates(microenvironment& m, diffusion_solver& d_solver, index_t index)
{
	if (!compute_internalized_substrates_)
		return;

	switch (m.mesh.dims)
	{
		case 1:
			release_dim<1>(d_solver.get_substrates_layout(), retrieve_agent_data(m.agents), m.mesh,
						   d_solver.get_substrates_pointer(), index);
			return;
		case 2:
			release_dim<2>(d_solver.get_substrates_layout(), retrieve_agent_data(m.agents), m.mesh,
						   d_solver.get_substrates_pointer(), index);
			return;
		case 3:
			release_dim<3>(d_solver.get_substrates_layout(), retrieve_agent_data(m.agents), m.mesh,
						   d_solver.get_substrates_pointer(), index);
			return;
		default:
			assert(false);
			return;
	}
}

void cell_solver::resize(const microenvironment& m)
{
	auto prev_capacity = numerators_.capacity();

	numerators_.resize(m.substrates_count * m.agents.size());
	denominators_.resize(m.substrates_count * m.agents.size());
	factors_.resize(m.substrates_count * m.agents.size());

	auto new_capacity = numerators_.capacity();

	if (new_capacity != prev_capacity)
	{
		reduced_numerators_ = std::make_unique<std::atomic<real_t>[]>(new_capacity);
		reduced_denominators_ = std::make_unique<std::atomic<real_t>[]>(new_capacity);
		reduced_factors_ = std::make_unique<std::atomic<real_t>[]>(new_capacity);
	}
}

void cell_solver::initialize(const microenvironment& m)
{
	compute_internalized_substrates_ = m.compute_internalized_substrates;

	resize(m);

	ballots_ = std::make_unique<std::atomic<index_t>[]>(m.mesh.voxel_count());
}
