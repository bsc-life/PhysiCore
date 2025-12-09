---
layout: default
title: Repository Structure
nav_order: 4
description: "Code organization guide with directory layout, module structure, and navigation tips"
---

# Repository Structure

This page documents the directory layout and organization of the PhysiCore repository. Understanding this structure will help you navigate the codebase and find the components you need.

## 1. Modules

PhysiCore is organized into modular libraries, each representing a different time-scale or aspect of multicellular simulation. All modules follow a consistent structure with clear separation between public APIs and implementation details.

### Module Directories

The following top-level directories contain PhysiCore modules:

- **`common/`** - Core abstractions and interfaces
- **`reactions-diffusion/`** - Diffusion and reaction solvers
- **`mechanics/`** - Mechanical interaction models
- **`phenotype/`** - Phenotype integration and orchestration

### Standard Module Structure

Each module follows this canonical structure:

```
module-name/
├── CMakeLists.txt         # Build configuration
├── include/               # Public API headers (exported via FILE_SET)
│   └── module-name/
│       ├── interface.h
│       └── types.h
├── src/                   # Private implementation files
│   ├── implementation.cpp
│   └── internal.h
├── tests/                 # Unit tests (GoogleTest)
│   ├── CMakeLists.txt
│   └── test_feature.cpp
├── kernels/               # Optional: hardware-specific backends
│   └── backend-name/
│       ├── CMakeLists.txt
│       ├── include/
│       └── src/
└── examples/              # Optional: demonstration applications
    ├── CMakeLists.txt
    └── example.cpp
```

**Key conventions:**
- **`include/`** - Public stable API, exported via CMake `FILE_SET HEADERS`
- **`src/`** - Private implementation, not visible to consumers
- **`tests/`** - Unit tests for the module
- **`kernels/`** - Self-contained solver backends (for pluggable implementations)
- **`examples/`** - Sample applications demonstrating usage

## 2. vcpkg - Dependency Management

PhysiCore uses [vcpkg](https://vcpkg.io/) for reproducible C++ dependency management. The following files and directories control dependency resolution:

### vcpkg Files and Directories

```
PhysiCore/
├── vcpkg/                       # vcpkg tool (git submodule)
├── vcpkg.json                   # Dependency manifest
├── vcpkg-configuration.json     # Baseline and overlay configuration
├── ports-overlays/              # Custom port definitions
│   ├── cccl/
│   │   ├── portfile.cmake
│   │   └── vcpkg.json
│   ├── noarr-structures/
│   │   ├── portfile.cmake
│   │   └── vcpkg.json
│   └── vtk-ioxml/
│       ├── portfile.cmake
│       └── vcpkg.json
└── triplets-overlays/           # Custom platform configurations
    ├── arm64-osx-gcc.cmake
    └── x64-osx-gcc.cmake
```

### File Descriptions

- **`vcpkg/`** - Git submodule pointing to Microsoft's [vcpkg repository](https://github.com/microsoft/vcpkg). Contains the vcpkg package manager tool itself. Bootstrap with `./vcpkg/bootstrap-vcpkg.sh` (Unix) or `vcpkg\bootstrap-vcpkg.bat` (Windows).

- **`vcpkg.json`** - **Manifest file** declaring all dependencies (e.g., `gtest`, `highway`, `tbb`, `fmt`, `nlohmann-json`). CMake reads this file and automatically installs dependencies during configuration via manifest mode.

- **`vcpkg-configuration.json`** - Configures the vcpkg **baseline** (specific commit/version of the vcpkg registry) and points to **overlay directories** for custom ports and triplets. Ensures reproducible dependency versions across all builds.

- **`ports-overlays/`** - Contains **custom or modified vcpkg port definitions**. Each subdirectory represents a package:
  - `portfile.cmake` - Build instructions for the package
  - `vcpkg.json` - Package metadata and dependencies

  Use overlays when you need a dependency not in the official registry or require custom build flags.

- **`triplets-overlays/`** - Defines **custom triplet files** (platform-specific build configurations). For example, `arm64-osx-gcc.cmake` configures GCC for macOS ARM64 instead of the default Clang. Triplets control compiler, linker, and toolchain settings.

### How It Works

1. CMake reads `vcpkg.json` and `vcpkg-configuration.json`
2. vcpkg installs dependencies to `build/<preset>/vcpkg_installed/`
3. Custom ports from `ports-overlays/` override official ports
4. Custom triplets from `triplets-overlays/` define build configurations
5. Dependencies are automatically linked via CMake's `find_package()`

## 3. Code Quality Tools

### `.clang-format`

Enforces consistent code formatting across the entire codebase.

**Purpose:** Defines formatting rules for C++ code (indentation, braces, spacing, line wrapping, etc.). PhysiCore uses a custom style based on industry best practices.

**Usage:**
```bash
# Format all files in the repository
git ls-files -z -- '*.h' '*.hpp' '*.cpp' '*.cu' | xargs -0 clang-format -i -style=file

# Or use the provided task
# VS Code: Run Task > Fix formatting with clang-format
```

**Integration:**
- GitHub Actions workflow `formal-checks.yml` validates formatting on every PR
- Pre-commit hooks can auto-format on commit (optional)
- IDE support: VS Code and CLion read `.clang-format` automatically

### `.clang-tidy`

Static analysis configuration for catching bugs, enforcing best practices, and ensuring code quality.

**Purpose:** Configures clang-tidy checks for:
- Modernization (use C++20 features)
- Bug detection (null pointer dereferences, memory leaks)
- Performance issues (unnecessary copies, inefficient algorithms)
- Readability and maintainability

**Usage:**
```bash
# Run linting on all files
run-clang-tidy -p build

# Auto-fix issues
run-clang-tidy -p build -fix

# Or use the provided tasks
# VS Code: Run Task > Run linting with clang-tidy
```

**Integration:**
- GitHub Actions workflow `formal-checks.yml` runs clang-tidy on every PR
- Requires a build directory with `compile_commands.json` (generated by CMake)
- IDE support: clangd language server uses `.clang-tidy` for diagnostics

## 4. Development Containers

### `.devcontainer/`

Defines containerized development environments for consistent, reproducible builds across all platforms.

**Purpose:** Provides a fully-configured Docker container with all required tools pre-installed:
- Compilers: GCC, Clang, MSVC (via Wine on Linux)
- Build tools: CMake, Ninja, vcpkg
- Code quality: clang-format, clang-tidy, clangd
- Debugging: gdb, lldb
- Sanitizers: ASAN, LSAN, UBSAN, TSAN

**Structure:**
```
.devcontainer/
├── devcontainer.json      # VS Code dev container configuration
├── Dockerfile             # Container image definition
└── setup.sh               # Post-create setup script
```

**Usage:**
1. Install [VS Code](https://code.visualstudio.com/) and [Docker](https://www.docker.com/)
2. Install the [Dev Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)
3. Open PhysiCore in VS Code
4. Click "Reopen in Container" when prompted

**Benefits:**
- Identical environment for all developers (eliminates "works on my machine")
- Pre-configured toolchain and dependencies
- Isolated from host system (no conflicts with system libraries)
- Fast onboarding for new contributors

## 5. GitHub Workflows

### `.github/workflows/`

Automated continuous integration and deployment pipelines.

```
.github/
├── workflows/
│   ├── cmake-multi-platform.yml      # Cross-platform build tests
│   ├── cmake-ubuntu-sanitized.yml    # Sanitizer builds (ASAN, TSAN, LSAN, UBSAN)
│   ├── formal-checks.yml             # Linting and formatting validation
│   ├── release.yml                   # Semantic versioning and releases
│   └── sonarqube.yml                 # Code quality analysis
└── actions/
    ├── cuda-toolkit/                 # Reusable CUDA setup action
    └── setup-vcpkg/                  # Reusable vcpkg caching action
```
