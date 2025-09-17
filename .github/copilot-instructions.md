<!--
These instructions are read by AI coding assistants to help them be productive in the
PhysiCore repository. Keep this file short and strictly actionable. Do not add
policy text or generic advice.
-->

# Copilot / AI assistant guidance — PhysiCore

Summary
- PhysiCore is a modular C++ simulation core. Major components (see top-level
  `CMakeLists.txt`) are:
  - `common/` — data structures, core interfaces (e.g. `process`, `base_agent`,
    `base_agent_data`). Look under `common/include` for canonical patterns.
  - `reactions-diffusion/biofvm/` — microenvironment / diffusion backends.
  - `mechanics/physicell/` — mechanics engines (implements `process` interface).
  - `phenotype/physicore/` — example executable that wires components together.

Big-picture/Why
- The code separates simulation core (data + interfaces) from pluggable
  numerical backends. Libraries are created per-subsystem and then linked into
  the example executable (`physicore.phenotype.physicore`). This lets the
  project swap diffusion / mechanics engines without changing core agent code.

Build, test, and debug (canonical commands)
- This repo uses CMake presets (see `CMakePresets.json`) and Ninja. Preferred
  workflows are via CMake presets. Example commands you can run in the repo root:

  - Configure + build (GCC Debug):
    cmake --preset=gcc-debug
    cmake --build --preset=gcc-debug

  - Configure + build + test (single preset):
    cmake --workflow=gcc-debug

  - Run tests for a preset (after configure/build):
    ctest --preset gcc-debug --output-on-failure

Notes about presets:
- Useful presets in `CMakePresets.json`: `gcc-debug`, `gcc-relwithdebinfo`,
  `gcc-release`, `llvm-debug`, `llvm-cov-debug`. The presets set sanitizer flags
  for Debug/TSAN/Coverage variants — prefer these when reproducing CI.

Project-specific conventions
- C++20, header-only public headers are exported using CMake FILE_SET HEADERS
  (see subproject CMakeLists). Use `physicore::...` target aliases when linking.
- Public API lives in `*/include/*.h`. Implementation belongs in `*/src/`.
- Tests are enabled with `PHYSICORE_BUILD_TESTS` (default on top-level builds).
  Tests use GoogleTest; tests live in each submodule under `tests/` and follow
  the pattern `*.tests` targets in CMake.
- Interfaces use small abstract classes in `common/include` (e.g., `process`).
  New engines/components should implement these interfaces and register as
  separate libraries that link `physicore::common`.

Integration points and examples
- Example wiring: `phenotype/physicore/CMakeLists.txt` and
  `phenotype/physicore/src/main.cpp` (executable entry) link
  `physicore::mechanics::physicell` and `physicore::reactions-diffusion::biofvm`.
- Example interface: `mechanics/physicell/include/environment.h` implements
  `process` (see `common/include/process.h`). Use `run_single_timestep()` as the
  per-step hook.

Quick tips for code edits
- When adding a new subsystem:
  1. Create a subdirectory and CMake target `add_library(physicore.<subsystem>)`.
  2. Export public headers with FILE_SET HEADERS under `include/` and add an
     ALIAS target `physicore::<subsystem>`.
  3. Link to `physicore::common` and follow the tests pattern if adding tests.
- Prefer modifying `common/include` for cross-cutting interfaces. Keep
  implementation details in `src/` to avoid exposing internal headers.

Files to read first (in order)
1. `README.md` — project intent and vision
2. `CMakePresets.json` — canonical build/test presets used by CI
3. `CMakeLists.txt` (top-level) — module boundaries
4. `common/include/*` — core interfaces and data shapes
5. `phenotype/physicore/src/main.cpp` — how components are wired at runtime

If something is ambiguous
- Ask the user to run a specific preset (`cmake --workflow=gcc-debug`) or to
  share failing build/test logs. For design questions, point to the small set
  of example files above so reviewers can see real usage.

End — ask for feedback if anything is missing or unclear.
