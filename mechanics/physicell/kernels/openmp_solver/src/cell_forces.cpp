#include "physicell/openmp_solver/cell_forces.h"

#include <algorithm>
#include <cmath>

namespace physicore::mechanics::physicell::kernels::openmp_solver {

namespace {

constexpr real_t kMinDistance = static_cast<real_t>(1e-5);
constexpr real_t kSimplePressureCoefficient = static_cast<real_t>(1.0);

} // namespace



} // namespace physicore::mechanics::physicell::kernels::openmp_solver
