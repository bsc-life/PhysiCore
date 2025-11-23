# Contributing to PhysiCore

Thank you for considering a contribution to PhysiCore. Your ideas, fixes, tests, and documentation improvements help shape a robust, modular C++20 simulation framework for multicellular biology. This guide explains what we are looking for, how to get started quickly, and how we review and accept changes.

## Table of Contents
- [Introduction](#introduction)
- [Ground Rules](#ground-rules)
- [Your First Contribution](#your-first-contribution)
- [Getting Started](#getting-started)
  - [Repository Layout Recap](#repository-layout-recap)
  - [Environment & Build](#environment--build)
  - [Making a Change (Standard Flow)](#making-a-change-standard-flow)
  - [Small or "Obvious" Fix Policy](#small-or-obvious-fix-policy)
  - [Code Style](#code-style)
  - [Commit Message Convention](#commit-message-convention)
  - [Adding a New Subsystem / Library](#adding-a-new-subsystem--library)
  - [Testing Guidelines](#testing-guidelines)
- [How to Report a Bug](#how-to-report-a-bug)
- [How to Suggest a Feature or Enhancement](#how-to-suggest-a-feature-or-enhancement)
- [Getting Help](#getting-help)
- [Acknowledgements](#acknowledgements)

## Introduction
PhysiCore is organized as small, composable libraries (mechanics, reactions-diffusion, phenotype, common) with stable public headers. Following these guidelines respects maintainers' time and increases the likelihood that your change will be merged promptly.

We welcome (and actively encourage) contributions across:
- New diffusion solver backends (BioFVM kernels)
- Mechanics features & environment extensions
- Performance improvements (SIMD, parallelism, memory layout)
- Testing (unit, integration, regression cases)
- Documentation & tutorial examples
- Developer tooling & CI improvements

We also value non-code contributions: issue triage, reproductions, clarifying docs, and architectural discussion.

## Ground Rules
Please read our [Code of Conduct](CODE_OF_CONDUCT.md) before contributing. Core expectations:
- Cross-platform consideration: avoid OS-specific assumptions; CI covers Linux, macOS, Windows.
- Tests: additions or changes to logic must include or extend tests. Prefer focused, deterministic tests.
- API stability: public headers under any `*/include/` folder form the published interfaces; coordinate breaking changes via an issue first.
- Modularity: keep new code in the appropriate subsystem; avoid unnecessary coupling.
- Performance: profile significant changes; do not prematurely micro-optimize without evidence.
- Inclusivity: be welcoming, patient, and constructive. Encourage first-time contributors.
- Small PRs: incremental, well-scoped changes are easier to review and merge.

## Your First Contribution
New here? Look for issues labeled `good first issue` and `help wanted`. These are curated to be approachable (small surface area, clear acceptance criteria). Unsure where to start? Open a discussion or issue describing what interests you and we can help scope something.

Resources for first-time open source contributors:
- https://makeapullrequest.com/
- https://www.firsttimersonly.com/
- Free video series: "How to Contribute to an Open Source Project on GitHub" (egghead.io)

Feel free to ask questions—everyone starts somewhere.

## Getting Started
### Repository Layout Recap
Repository is divided into directories according to the typical time-scale of agent-based simulator.

- `common/`: core abstractions (`timestep_executor.h`, shared types, utilities)
- `reactions-diffusion/`: diffusion solvers.
- `mechanics/`: mechanics interaction between agents.
- `phenotype/`: phenotypicall properties of cells.

### Environment & Build
Requirements: C++20 compiler (GCC ≥ 13 / Clang ≥ 17 / MSVC recent), CMake ≥ 3.27 with preset support, Ninja, Git.

Bootstrapping vcpkg (recommended local instance):
```sh
export VCPKG_ROOT="$PWD/vcpkg"
"$VCPKG_ROOT"/bootstrap-vcpkg.sh
```
Configure, build, test (example GCC debug flow):
```sh
cmake --preset=gcc-debug
cmake --build --preset=gcc-debug
ctest --preset gcc-debug --output-on-failure
```
Other presets (LLVM, Release, sanitizer builds) live in `CMakePresets.json`. Use `cmake --workflow=<preset>` for a one-shot configure/build/test.

Dev Container: Opening the repo in VS Code with the dev container yields a pre-installed toolchain and presets. This is the fastest path to a green build. We currently maintain two dev containers: `physicore-dev` and `physicore-cuda-dev` installing general and CUDA related dependencies respectively.

### Making a Change (Standard Flow)
1. Open / find or create an issue describing the problem/feature; gather consensus for larger changes.
2. Create a feature branch (`git checkout -b feature/short-slug`).
3. Implement incrementally; keep commits logically grouped.
4. Add / update tests (new behavior must be tested; bug fixes should include regression tests).
5. Run the full debug test preset locally; optionally run sanitizer builds for memory/thread safety.
6. Ensure formatting (clang-format) and no warnings you can reasonably fix.
7. Open a Pull Request (PR) referencing the issue. Fill in the PR template completely.

### Small or "Obvious" Fix Policy
Typo corrections, minor spelling, whitespace, and purely internal comment adjustments can be submitted without prior issue discussion. These must not alter public behavior or API. If unsure, open an issue first.

### Code Style
- C++20, favor standard library & `<chrono>`, `<span>`, `<ranges>` where appropriate.
- RAII for resource management; avoid raw `new` / `delete` (use smart pointers or value types).
- Prefer `enum class` over plain `enum`.
- Avoid macros; use `constexpr` or templates.
- Keep functions short & cohesive; refactor large logic blocks.
- Use GoogleTest for unit/integration tests; follow existing test patterns.
- Format code with `clang-format` (see `.clang-format`).

### Commit Message Convention
We use [Semantic Versioning](https://semver.org/) driven by [Conventional Commits](https://www.conventionalcommits.org/). Commit messages MUST follow the following format. Our automated release tooling relies on this structure to generate changelogs and determine version bumps.
```
<type>(scope): summary

[optional body]

[optional footer(s)]
```

Examples:
- `feat(biofvm): add OpenMP diffusion kernel`
- `fix(physicell): correct boundary force accumulation`

### Adding a New Subsystem / Library
Follow checklist (see `copilot-instructions.md` for detailed steps):
1. Create `subsystem/{include,src}` and `CMakeLists.txt` with `add_library` + alias `physicore::subsystem`.
2. Export headers via `FILE_SET HEADERS`.
3. Link against `physicore::common` (and others as needed).
4. Add tests under `subsystem/tests` using GoogleTest.
5. Document purpose & integration in `README` or subsystem-specific docs.

### Testing Guidelines
- Unit tests: fast, deterministic; cover edge cases & error paths.
- Integration tests: validate solver-mechanics interplay; keep runtime bounded.
- Benchmarks (if added): isolated, optional (may be excluded from default test preset).
- Run: `ctest --preset gcc-debug --output-on-failure` (adapt preset as needed).

## How to Report a Bug
Open a GitHub issue including:
1. PhysiCore version or commit SHA (`git rev-parse HEAD`).
2. Toolchain details (compiler & version, OS, architecture).
3. Exact steps or minimal code snippet to reproduce.
4. Expected behavior.
5. Actual behavior (include logs, error messages, stack trace if applicable).
6. Any performance implications (e.g., slowdown factor, memory growth).

For runtime crashes, include backtrace (`gdb` or `lldb`), and compiler flags used. For diffusion/mechanics numerical issues, share parameter set and initial conditions.

## How to Suggest a Feature or Enhancement
Open an issue labeled `enhancement` describing:
- Motivation & problem statement (why current behavior insufficient)
- Proposed solution (API surface, data structures, algorithms)
- Alternatives considered
- Potential performance / memory impact
- Testing strategy (unit vs integration benchmarks)

Architectural philosophy: keep subsystems independent with clearly defined interfaces; prefer opt-in for heavyweight dependencies; maximize testability. New solver backends should register via the BioFVM solver registry pattern.

## Getting Help
If stuck, open a new [discussion](https://github.com/bsc-life/PhysiCore/discussions) and describe your issue. Community members and maintainers monitor discussions and can provide guidance.

## Acknowledgements
Your effort—large or small—matters. Thank you for helping advance PhysiCore.

---
Happy hacking!
