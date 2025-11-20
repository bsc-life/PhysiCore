#pragma once

// Forward declarations for both TBB and CUDA thrust solver backends
namespace physicore::biofvm::kernels::tbb_thrust_solver {
void attach_to_registry();
}

namespace physicore::biofvm::kernels::cuda_thrust_solver {
void attach_to_registry();
}
