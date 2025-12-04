#pragma once

#include <array>
#include <unordered_map>
#include <vector>

#include "spatial_index.h"

namespace physicore::mechanics::micromechanics {

struct grid_key
{
	int x, y, z;
	bool operator==(const grid_key& other) const { return x == other.x && y == other.y && z == other.z; }
};

struct grid_key_hasher
{
	std::size_t operator()(const grid_key& k) const
	{
		// Simple hash combining
		return ((std::hash<int>()(k.x) ^ (std::hash<int>()(k.y) << 1)) >> 1) ^ (std::hash<int>()(k.z) << 1);
	}
};

class uniform_grid_spatial_index : public spatial_index
{
	real_t cell_size;
	std::unordered_map<grid_key, std::vector<index_t>, grid_key_hasher> grid;

public:
	explicit uniform_grid_spatial_index(real_t cell_size = 30.0); // Default to typical cell size

	void build(const environment& env) override;
	std::vector<index_t> query_neighbors(const environment& env, index_t agent_index, real_t radius) const override;
};

} // namespace physicore::mechanics::micromechanics
