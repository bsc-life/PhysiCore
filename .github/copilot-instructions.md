<!--
These instructions are read by AI coding assistants to help them be productive in the
PhysiCore repository. Keep this file short and strictly actionable. Do not add
policy text or generic advice.
-->

# Copilot / AI assistant guidance — PhysiCore

Summary
- Modular C++20 simulation core. Top-level `CMakeLists.txt` adds: `common/`, `reactions-diffusion/biofvm/`, `mechanics/physicell/`, `phenotype/physicore/`.
- Stable core API in `common/include` with minimal interfaces; numerical engines live as independent libraries and are wired into an example executable.

Architecture essentials (what talks to what)
- Core timestep API: `common/include/process.h` defines `class timestep_executor { run_single_timestep(); serialize_state(); }`.
- Mechanics implements the timestep executor: `mechanics/physicell/include/environment.h` with `environment::run_single_timestep()` and `serialize_state()` in `src/environment.cpp`.
- Microenvironment (BioFVM) provides diffusion, serializers, and solver plug-ins:
  - Serializer API: `reactions-diffusion/biofvm/include/serializer.h`; VTK serializer in `include/vtk_serializer*.h` + `src/vtk_serializer*.cpp`. It writes `.vti` files and a collection `.pvd` under the chosen output dir (tests use `vtk_microenvironment/`).
  - Solver registry for pluggable backends: `include/solver_registry.h` with `registry_adder` pattern. Kernels register via `kernels/*/src/register_solver.cpp`; all kernels attached in `src/solver_registry_attach.cpp` (gated by `PHYSICORE_HAS_THRUST`).
- Wiring: the example app `phenotype/physicore` links `physicore::mechanics::physicell` and `physicore::reactions-diffusion::biofvm`.

Build, test, debug (use presets exactly like CI)
- Configure + build (GCC debug):
  - cmake --preset=gcc-debug
  - cmake --build --preset=gcc-debug
- Full cycle: cmake --workflow=gcc-debug
- Run tests: ctest --preset gcc-debug --output-on-failure
- Other presets: llvm-{debug,release,relwithdebinfo}, gcc-{release,relwithdebinfo}, appleclang-*, msvc-* (see `CMakePresets.json`). Sanitizers and TSAN suppressions come from repo root (`asan-suppressions.txt`, `tsan-suppressions.txt`).
- vcpkg: CMake uses `$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake` and the manifest `vcpkg.json` (includes VTK IOXML).

Project conventions
- Public headers live in `*/include/` and are exported via CMake `FILE_SET HEADERS`; link against aliases like `physicore::mechanics::physicell` or `physicore::reactions-diffusion::biofvm`.
- Tests are on by default for top-level builds via `PHYSICORE_BUILD_TESTS`; per-subproject tests live under `*/tests/` and use GoogleTest with discovered `*.tests` executables.
- Kernels live in `*/kernels/<backend>/` and expose an `attach_to_registry()` entrypoint; `src/solver_registry_attach.cpp` ensures they self-register at link/load.

Concrete pointers and examples
- Implementing the timestep hook: see `mechanics/physicell/include/environment.h` and `src/environment.cpp` for `run_single_timestep()`.
- Adding a solver backend: follow `reactions-diffusion/biofvm/kernels/openmp_solver/src/register_solver.cpp` and the `registry_adder` usage in `include/solver_registry.h`.
- Writing VTK output: `reactions-diffusion/biofvm/include/vtk_serializer.h` + `src/vtk_serializer.cpp`; base class handles `.pvd` aggregation in `src/vtk_serializer_base.cpp`.
- Wiring an app: `phenotype/physicore/CMakeLists.txt` links mechanics and RD libraries; `src/main.cpp` is the integration point.

Quick checklist to add a subsystem/library
1) Create `your_subsystem/{include,src}` and add `add_library(physicore.your_subsystem)`.
2) Export headers via `FILE_SET HEADERS`; add alias `physicore::your_subsystem`.
3) Link to `physicore::common` and add `tests/` with a `*.tests` target using GTest.
4) If runtime-pluggable, mirror the registry pattern shown in BioFVM.

Files to skim first
- `README.md`, `CMakePresets.json`, top-level `CMakeLists.txt`, `common/include/*`, `phenotype/physicore/src/main.cpp`.

If something is ambiguous
- Provide the preset and exact build/test logs (e.g., `cmake --workflow=gcc-debug`); most wiring lives in the files listed above.

End — please suggest any missing conventions or workflows you rely on so we can capture them here.
