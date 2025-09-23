<!--
These instructions are read by AI coding assistants to help them be productive in the
PhysiCore repository. Keep this file short and strictly actionable. Do not add
policy text or generic advice.
-->

# Copilot / AI assistant guidance — PhysiCore

Summary
- PhysiCore is a modular C++20 simulation core. Major components (see top-level
  `CMakeLists.txt`):
  - `common/` — core data shapes and small abstract interfaces (look in
    `common/include`). Example: `timestep_executor` with `run_single_timestep()`.
  - `reactions-diffusion/biofvm/` — microenvironment and diffusion solver code.
  - `mechanics/physicell/` — mechanics engine implementing the timestep API.
  - `phenotype/physicore/` — example executable wiring subsystems together.

Big-picture (why it looks like this)
- The repo keeps a tiny, stable core API in `common/` and implements numerical
  engines as independently-buildable libraries. Subsystems are linked into the
  example executable so implementations can be swapped without changing core
  agent logic.

Build / test / debug (exact commands)
- This project uses CMake presets + Ninja. Use presets to reproduce CI precisely.
  - Configure + build (GCC debug):
    cmake --preset=gcc-debug
    cmake --build --preset=gcc-debug
  - Full configure + build + test workflow:
    cmake --workflow=gcc-debug
  - Run tests (after configure/build):
    ctest --preset gcc-debug --output-on-failure
- Notes: presets configure sanitizers and toolchain via `$VCPKG_ROOT`. Ensure
  vcpkg is available or `VCPKG_ROOT` is set. TSAN builds reference
  `tsan-suppressions.txt` from the repo root.

Project conventions and patterns
- Language: C++20. Public API: `*/include/*.h`. Implementation: `*/src/`.
- Subprojects export headers using CMake `FILE_SET HEADERS` and expose alias
  targets `physicore::<subsystem>` (see subproject CMakeLists in
  `mechanics/physicell/` and `reactions-diffusion/biofvm/`). Link to those
  aliases from `phenotype/physicore` when wiring an executable.
- Tests: enabled by `PHYSICORE_BUILD_TESTS` (defaults ON for top-level builds).
  Tests live under `*/tests/` and use GoogleTest with CMake `*.tests` targets.
- Interfaces: small abstract classes live in `common/include`. Example:
  `common/include/process.h` defines `timestep_executor::run_single_timestep()`;
  implementors in `mechanics/physicell/` or `reactions-diffusion/biofvm/` should
  implement that method.

Integration points and concrete examples
- Wiring example: `phenotype/physicore/CMakeLists.txt` and
  `phenotype/physicore/src/main.cpp` show how libraries are linked into the
  executable — use them to see which subsystems are active by default.
- Interface example: `mechanics/physicell/include/environment.h` implements the
  core timestep interface. Use `run_single_timestep()` as the per-step hook.
- VCPKG: project CMake uses `$env{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake`.
  If CI fails due to missing packages, ensure vcpkg dependencies are installed
  via the repository `vcpkg.json` manifest.

Quick checklist for adding a subsystem
1. Create `your_subsystem/` with `include/` and `src/`.
2. Add `add_library(physicore.your_subsystem)` in that subproject CMakeLists.
3. Export headers with `FILE_SET HEADERS` and add an alias `physicore::your_subsystem`.
4. Link to `physicore::common` and add tests under `tests/` following `*.tests`.

Files to read first
1. `README.md` — high-level intent
2. `CMakePresets.json` — canonical build/test presets
3. `CMakeLists.txt` (top-level) — module boundaries
4. `common/include/*` — core interfaces and data shapes
5. `phenotype/physicore/src/main.cpp` — wiring example

If something is ambiguous
- Ask for the failing preset + build/test logs (run `cmake --workflow=gcc-debug`)
  or point reviewers to the small set of wiring files above.

End — request feedback on any missing or unclear sections so I can iterate.
