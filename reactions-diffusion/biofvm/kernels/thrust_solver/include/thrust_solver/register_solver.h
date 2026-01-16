#pragma once

// Forward declarations for both TBB and CUDA thrust solver backends
namespace physicore::reactions_diffusion::biofvm::kernels::tbb_thrust_solver {
void attach_to_registry();
}

namespace physicore::reactions_diffusion::biofvm::kernels::cuda_thrust_solver {
void attach_to_registry();
}
