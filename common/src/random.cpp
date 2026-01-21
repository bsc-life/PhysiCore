#include <omp.h>
#include <random>
#include <vector>

#include <common/random.h>


thread_local std::mt19937 generator;

physicore::random& physicore::random::instance()
{
	static random instance;
	return instance;
}

physicore::real_t physicore::random::uniform(const physicore::real_t min, const physicore::real_t max)
{
	std::uniform_real_distribution<physicore::real_t> distribution(min, max);
	return distribution(generator);
}

physicore::real_t physicore::random::normal(const physicore::real_t mean, const physicore::real_t std)
{
	std::normal_distribution<physicore::real_t> distribution(mean, std);
	return distribution(generator);
}

void physicore::random::set_seed(unsigned int seed)
{
#ifdef _OPENMP
	std::vector<unsigned int> initial_sequence(omp_get_num_threads());

	for (int i = 0; i < omp_get_num_threads(); i++)
		initial_sequence[i] = seed + i;

	std::seed_seq seq(initial_sequence.begin(), initial_sequence.end());

	std::vector<unsigned int> seeds(omp_get_num_threads());
	seq.generate(seeds.begin(), seeds.end());

	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		generator.seed(seeds[id]);
	}
#else
	generator.seed(seed);
#endif
}
