# PhysiCore AI Agent Instructions

PhysiCore is a modular C++20 framework for agent-based multicellular simulations, re-architecting PhysiCell with pluggable diffusion solvers and mechanics engines.

## Architecture Overview

**Modular Design**: The codebase is divided into time-scale-based subsystems, each a separate library with stable public interfaces via CMake `FILE_SET HEADERS`:

- `common/` - Core abstractions (`timestep_executor`, agent interfaces, types)
- `reactions-diffusion/biofvm/` - Diffusion/reaction solvers with self-registering backends
- `mechanics/physicell/` - Mechanics engines and force models
- `phenotype/physicore/` - Phenotype models and executable wiring

**Public vs Internal APIs**: Headers in `*/include/` are public stable APIs. Implementation details in `*/src/` are private. All libraries provide namespace aliases (e.g., `physicore::common`, `physicore::reactions-diffusion::biofvm`).

**Pluggable Solver Pattern**: The `solver_registry` enables runtime selection of diffusion backends:
- Backends self-register via `registry_adder<SolverT>` template in their registration units
- See `reactions-diffusion/biofvm/kernels/{openmp_solver,thrust_solver}/src/register_solver.cpp`
- Registration happens at static initialization in `solver_registry_attach.cpp`
- Use `solver_registry::instance().get("solver_name")` to instantiate solvers

**Key Abstraction**: `timestep_executor` (in `common/timestep_executor.h`) defines the simulation loop contract:
```cpp
virtual void run_single_timestep() = 0;
virtual void serialize_state(real_t current_time) = 0;
```

## Build System (CMake + vcpkg)

**Presets-First Workflow**: Use CMake presets exclusively (defined in `CMakePresets.json`):
```sh
cmake --preset=gcc-debug          # Configure
cmake --build --preset=gcc-debug  # Build
ctest --preset gcc-debug          # Test
cmake --workflow --preset=gcc-debug  # All-in-one
```

Available presets: `{gcc,llvm,appleclang,msvc}-{debug,release,relwithdebinfo}`

**vcpkg Integration**: Dependencies managed via `vcpkg.json` manifest mode. The `VCPKG_ROOT` environment variable must point to vcpkg (repo includes it as submodule at `./vcpkg`). First configure automatically installs dependencies under `build/<preset>/vcpkg_installed/`.

**Adding Libraries**: Follow the CMake pattern in existing subsystems:
1. `add_library(physicore.subsystem.name)` + `add_library(physicore::subsystem::name ALIAS ...)`
2. Export public headers via `FILE_SET HEADERS BASE_DIRS include/`
3. Link `physicore::common` and other dependencies
4. Tests go in `subsystem/tests/` with GoogleTest

**Sanitizer Builds**: Debug presets include ASAN/LSAN/UBSAN; RelWithDebInfo includes TSAN. Suppressions in `lsan.supp`/`tsan.supp`.

## Development Workflows

**Testing**: Use GoogleTest with `TEST()` or `TEST_F()` macros. Run via `ctest --preset <preset> --output-on-failure`. Tests live alongside implementation in `*/tests/`.

**Commit Messages**: Strictly follow Conventional Commits (`<type>(scope): summary`). Examples: `feat(biofvm): add CUDA kernel`, `fix(physicell): boundary force bug`. This drives semantic-release versioning.

## Code Conventions

**C++20 Modern Style**:
- Use `std::span`, `std::ranges`, `<chrono>`, `constexpr` over macros
- RAII for resources (smart pointers, value types) - no raw `new`/`delete`
- `enum class` over plain `enum`
- Namespaces: `physicore` for common, `physicore::biofvm` for reactions-diffusion, nested `::kernels::` for backends

**Agent Data Pattern**: Agents use structure-of-arrays (SoA) via `base_agent_data` and proxy objects (`base_agent`) for element access. See `common/include/common/base_agent.h` and tests for usage.
