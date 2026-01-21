#include "common_solver.h"

#include <span>

namespace physicore::mechanics::physicell::kernels::openmp_solver {

std::size_t common_solver::get_mesh_index(const voxel_pos_t& position, const cartesian_mesh& mesh)
{
	return mesh.linearize(position[0], position[1], position[2]);
}

common_solver::voxel_pos_t common_solver::get_mesh_position(const real_t* position, const cartesian_mesh& mesh)
{
	return mesh.voxel_position(std::span<const real_t>(position, static_cast<std::size_t>(mesh.dims)));
}

} // namespace physicore::mechanics::physicell::kernels::openmp_solver
