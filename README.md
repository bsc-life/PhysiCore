[![CMake Build on Ubuntu, Windows, and MacOS](https://github.com/bsc-life/PhysiCore/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/bsc-life/PhysiCore/actions/workflows/cmake-multi-platform.yml)
[![ASAN, LSAN, UBSAN, TSAN Build on Ubuntu](https://github.com/bsc-life/PhysiCore/actions/workflows/cmake-ubuntu-sanitized.yml/badge.svg)](https://github.com/bsc-life/PhysiCore/actions/workflows/cmake-ubuntu-sanitized.yml)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=bsc-life_PhysiCore&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=bsc-life_PhysiCore)
[![Release](https://img.shields.io/github/v/release/bsc-life/PhysiCore)](https://github.com/bsc-life/PhysiCore/releases/latest)
[![semantic-release: angular](https://img.shields.io/badge/semantic--release-angular-e10079?logo=semantic-release)](https://github.com/semantic-release/semantic-release)

# PhysiCore

**PhysiCore** is a modern re-architecture of the [PhysiCell](http://physicell.org) agent-based simulation framework, built for flexibility, performance, and clarity.
It introduces a cleaner, modular codebase that supports **interchangeable diffusion and mechanics models**, making it easy to experiment with new numerical methods, physical formulations, and biological processes without overhauling the entire simulation engine.

---

## Key Features

- **Modernized core architecture** – streamlined, maintainable, and extensible.
- **Pluggable diffusion solvers** – swap between FEM, FVM, FDM, or custom schemes.
- **Interchangeable mechanics engines** – easily test alternative force models and integration methods.
- **Performance-friendly design** – ready for HPC, vectorization, and GPU offload.
- **Clear separation of concerns** – clean interfaces between core simulation logic and numerical backends.

---

## Vision

PhysiCore aims to be the go-to foundation for **next-generation, high-performance multicellular simulations**, enabling researchers to focus on science instead of wrestling with legacy code.

---

## Code Structure

```
PhysiCore/
|-- common/                  (Core interfaces; the timestep contract lives in timestep_executor.h)
|   \-- include/
|       \-- timestep_executor.h
|-- reactions-diffusion/     (Diffusion and reaction solvers with pluggable backends)
|   \-- biofvm/
|       |-- include/
|       |-- src/
|       |-- examples/
|       \-- kernels/         (Solver backends that self-register)
|-- mechanics/               (Mechanics engines and force models)
|   \-- ...
|-- phenotype/               (Phenotype models and wiring into executables)
|   \-- ...
|-- ports-overlays/          (vcpkg overlays for pinned third-party ports)
\-- build/                   (Preset-specific build trees with vcpkg_installed/)
```

### Typical Subproject Layout

- `include/` – Public headers exported via CMake `FILE_SET` entries so downstream targets can consume stable interfaces.
- `src/` – Library implementation files; often link against `physicore::common` and register concrete components.
- `examples/` – Optional sample applications demonstrating how to integrate the library pieces in isolation.
- `kernels/` – Solver backends or hardware-specific implementations that self-register through the subsystem's registry entry points.

## Dependencies and vcpkg

This project uses vcpkg in manifest mode. CMake presets automatically point to the vcpkg toolchain via the environment variable `VCPKG_ROOT`.

- If you already have vcpkg installed, ensure `VCPKG_ROOT` points to it.
- If not, you can use the repo-local copy under `./vcpkg` (recommended for reproducibility).

Set up on Linux/macOS (zsh):

```sh
# Option A: use the repo-local vcpkg
git submodule update --init
export VCPKG_ROOT="$PWD/vcpkg"
"$VCPKG_ROOT"/bootstrap-vcpkg.sh

# Option B: point to your existing vcpkg install
export VCPKG_ROOT="/path/to/your/vcpkg"
```

Set up on Windows (PowerShell):

```powershell
# Option A: repo-local vcpkg
git submodule update --init
$env:VCPKG_ROOT = "$PWD/vcpkg"
& "$env:VCPKG_ROOT/bootstrap-vcpkg.bat"

# Option B: existing vcpkg install
$env:VCPKG_ROOT = "C:\\path\\to\\vcpkg"
```

Manifest dependencies (from `vcpkg.json`):

- `cccl` - NVIDIA Thrust library for TBB and CUDA backends
- `gtest` - Testing framework
- `highway` - Vendor-agnostic vectorization library
- `noarr-structures` - Memory layouts library
- `pugixml` - Lightweight XML parsing library for configuration files
- `tbb` - Intel OneAPI Threading Building Blocks
- `vtk-ioxml` - VTK IO for BioFVM serializer

Notes
- The repo ships `vcpkg-configuration.json` with a pinned baseline and `ports-overlays/` for custom/overlay ports.
- On first configure, vcpkg will build and install dependencies for your preset under `build/<preset>/vcpkg_installed/`.


## Build

Prerequisites
- CMake with presets support, a C++20 compiler (GCC/Clang/MSVC), and Ninja (presets use the Ninja generator).
- vcpkg available via `VCPKG_ROOT` as described above.

Common flows (Linux/macOS, zsh):

```sh
# Configure + build (GCC release)
cmake --preset=gcc-release
cmake --build --preset=gcc-release

# Run tests
ctest --preset gcc-release --output-on-failure

# Full configure+build+test workflow in one step
cmake --workflow --preset=gcc-release
```

Other useful presets (see `CMakePresets.json`):
- `llvm-{debug,release}` for clang compilation
- `gcc-{debug,release}` for gcc compilation
- `appleclang-{debug,release}` for AppleClang compilation
- `msvc-{debug,release}` for Microsoft cl compilation

## Example Projects

The repository provides several example projects under `examples/` to demonstrate PhysiCore's capabilities:

### Diffussion

A simple example showcasing diffusion of a single substrate in a 3D domain with static cells acting as sources and sinks. It runs a 30m simulation and outputs VTK-compatible file.
To build and run:

```sh
cmake --preset=gcc-release
cmake --build --preset=gcc-release --target physicore.reactions-diffusion.biofvm.diffuse
build/gcc-release/reactions-diffusion/biofvm/examples/physicore.reactions-diffusion.biofvm.diffuse
```

The run generates `output/` with VTK files for visualization.

### Mechanics

The `mechanics/micromechanics` module provides cell mechanics: neighbor detection, force calculation (repulsion, adhesion, spring-damper), and position integration.

**Key Components:**

| Component | Role |
|-----------|------|
| `environment` | Central simulation state: agents, solver, spatial index, serializer |
| `agent_data` | SoA storage for agent properties (positions, velocities, radii, forces) |
| `cell_data` | Aggregated cell-level properties (computed from grouped agents) |
| `solver` | Force computation + position update (pluggable via `solver_registry`) |
| `simulation_parameters` | Type-based interaction rules + feature flags |

**Potentials:** `standard` (repulsion + adhesion), `morse` (smooth wells), `kelvin_voigt` (spring-damper).

**Minimal Example:**

```cpp
#include <micromechanics/environment.h>
#include <openmp_solver/register_solver.h>  // Registers "openmp_solver"

using namespace physicore::mechanics::micromechanics;

int main() {
    // Attach solver to registry (required for static library linking)
    kernels::openmp_solver::attach_to_registry();

    // Create environment with 0.01 timestep
    environment env(0.01);
    env.params.solver_name = "openmp_solver";
    env.initialize_solver();

    // Add agents via create() + configure properties
    auto* agent = env.agents->create();
    agent->position() = { 0.0, 0.0, 0.0 };
    agent->radius() = 10.0;

    // Simulation loop
    real_t t = 0.0, t_max = 10.0;
    while (t < t_max) {
        env.run_single_timestep();
        env.serialize_state(t);  // writes VTK to output/
        t += env.timestep;
    }
}
```

**Serialization:** VTK output is written to `output/` by default. The serializer creates `.vtu` files per timestep and a `.pvd` collection file for ParaView playback. Disable with `env.serializer.reset();`.
