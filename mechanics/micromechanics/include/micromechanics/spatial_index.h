#pragma once

#include <vector>

#include <common/types.h>

namespace physicore::mechanics::micromechanics {

class environment;

class spatial_index
{
public:
	virtual ~spatial_index() = default;

	virtual void build(const environment& env) = 0;
	virtual std::vector<index_t> query_neighbors(const environment& env, index_t agent_index, real_t radius) const = 0;
};

} // namespace physicore::mechanics::micromechanics
