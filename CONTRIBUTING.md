# Contributing to PhysiCore

Thanks for helping us build a modern, modular simulation framework. This guide outlines the expectations we have for contributors and the quickest way to get productive in the repository.

## Before You Start
- Read `README.md` to understand the repository layout, build presets, and dependency management with vcpkg.
- Familiarize yourself with the core subsystems (common, mechanics, reactions-diffusion, phenotype) so you know where a change should live.
- Search the existing issues and pull requests to avoid duplicating ongoing work.

## Environment Setup
- Ensure you have a C++20 compiler, CMake (with preset support), and Ninja.
- Point `VCPKG_ROOT` at the provided `./vcpkg` instance (recommended) or an existing vcpkg install, then bootstrap it:
  ```sh
  export VCPKG_ROOT="$PWD/vcpkg"
  "$VCPKG_ROOT"/bootstrap-vcpkg.sh
  ```
- Configure and build using the supplied presets. A common flow during development is:
  ```sh
  cmake --preset=gcc-debug
  cmake --build --preset=gcc-debug
  ctest --preset gcc-debug --output-on-failure
  ```
  Other presets such as `llvm-debug` or sanitizer-enabled builds are available in `CMakePresets.json`.

### VSCode Dev Container
Prefer the provided VS Code dev container for a pre-configured environment; it ships with the required compilers, CMake presets, vcpkg bootstrap, and clang-format so you can build and run tests immediately after opening the folder in VS Code.

## Start Contributing
Start by creating a feature branch from the default branch and keep each change focused. This branch will be the basis for your pull request. The contents of your pull request should be a logical unit that can be reviewed and tested independently, typically starting by defining a clear problem statement or issue. Before starting work, get familiar with the code structure and identify where your changes will fit.

### Coding Guidelines
- Follow the existing style of the subsystem you are touching (this project favors clear, modern C++20 with standard library facilities over custom utilities).
- Keep public interfaces stable: headers under `*/include/` constitute the exported API. Coordinate interface changes with maintainers.
- Maintain subsystem boundaries:
  - `common/` exposes core concepts and the timestep contract in `timestep_executor.h`.
  - `reactions-diffusion/biofvm/` owns diffusion kernels, serializers, and registry wiring.
  - `mechanics/physicell/` implements mechanics-specific timestep executors and environment logic.
  - `phenotype/physicore/` links the pieces into runnable examples or apps.
- Add succinct, informative comments only when the intent is not obvious from the code.

### Contributing Acceptance
To get your changes merged, please open a pull request against the default branch. Pay attention to the following:

- The pull request will be automatically built and tested on multiple platforms using GitHub Actions. Be prepared that your changes are supposed to work across Linux, Windows, and macOS.
- Please format your C++ sources with clang-format to maintain a consistent style. The repository includes a `.clang-format` file at the root. You can run clang-format manually or set up your editor to format on save.
- Your changes will be analyzed by SonarCloud for static analysis. Address any issues raised or provide justification if a finding is a false positive.

---

We appreciate every contribution, from bug fixes and documentation improvements to new solvers and mechanics models. Thank you for investing your time in PhysiCore!
