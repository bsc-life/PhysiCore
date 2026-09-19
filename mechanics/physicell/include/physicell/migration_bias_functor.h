#pragma once

#include <memory>

#include <common/types.h>

namespace physicore::mechanics::physicell {

enum class migration_bias_type
{
	none,
	simple,
	advanced,
	advanced_normalized
};

class migration_bias_functor
{
public:
	migration_bias_functor() = default;
	migration_bias_functor(const migration_bias_functor&) = delete;
	migration_bias_functor(migration_bias_functor&&) = delete;
	migration_bias_functor& operator=(const migration_bias_functor&) = delete;
	migration_bias_functor& operator=(migration_bias_functor&&) = delete;

	virtual void update_migration_bias(index_t agent_index, real_t* agent_migration_bias) = 0;
	virtual ~migration_bias_functor() = default;
};

using migration_bias_func_ptr = std::unique_ptr<migration_bias_functor>;

} // namespace physicore::mechanics::physicell
