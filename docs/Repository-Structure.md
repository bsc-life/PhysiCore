---
layout: default
title: Repository Structure
nav_order: 4
description: "Code organization guide with directory layout, module structure, and navigation tips"
---

# Repository Structure

This page documents the directory layout and organization of the PhysiCore repository. Understanding this structure will help you navigate the codebase and find the components you need.

## Top-Level Directory Layout

```
PhysiCore/
├── .devcontainer/         # Development container configuration
├── .github/               # GitHub Actions workflows and issue templates
├── .vscode/               # VS Code workspace settings
├── common/                # Core abstractions and interfaces
├── reactions-diffusion/   # Diffusion and reaction modules
├── mechanics/             # Mechanical interaction modules  
├── phenotype/             # Phenotype integration modules
├── ports-overlays/        # Custom vcpkg port overlays
├── triplets-overlays/     # Custom vcpkg triplet definitions
├── vcpkg/                 # vcpkg submodule (package manager)
├── build/                 # Build output (generated, not in git)
├── CMakeLists.txt         # Root CMake configuration
├── CMakePresets.json      # CMake preset definitions
├── vcpkg.json             # vcpkg manifest (dependencies)
├── vcpkg-configuration.json  # vcpkg baseline and overlays
├── .clang-format          # Code formatting rules
├── .clang-tidy            # Static analysis configuration
├── lsan.supp              # Leak sanitizer suppressions
├── tsan.supp              # Thread sanitizer suppressions
├── README.md              # Project overview
├── CONTRIBUTING.md        # Contribution guidelines
├── CODE_OF_CONDUCT.md     # Community code of conduct
└── LICENSE                # License information
```

## Module Structure

Each module (common, reactions-diffusion, mechanics, phenotype) follows a consistent structure. Here we describe the typical layout using examples from the repository.

### Common Module

The `common` module provides foundational types and interfaces used throughout PhysiCore.

```
common/
├── CMakeLists.txt         # Build configuration for common library
├── include/               # Public API headers
│   └── common/
│       ├── timestep_executor.h        # Main simulation loop contract
│       ├── base_agent.h               # Agent proxy interface
│       ├── base_agent_data.h          # SoA data structures
│       ├── base_agent_container.h     # Agent collection
│       ├── base_agent_interface.h     # Agent interface definition
│       ├── generic_agent_container.h  # Generic container
│       ├── generic_agent_solver.h     # Generic solver
│       ├── concepts.h                 # C++20 concepts
│       └── types.h                    # Common type definitions
└── tests/                 # Unit tests
    ├── CMakeLists.txt
    ├── base_agent_test.cpp
    └── ...
```

**Key Points:**
- Headers in `include/common/` are the **public API**
- Exported via CMake `FILE_SET HEADERS`
- All other modules depend on `physicore::common`

### Reactions-Diffusion Module

The reactions-diffusion module contains implementations for substrate transport and reaction kinetics.

```
reactions-diffusion/
└── biofvm/                # BioFVM implementation
    ├── CMakeLists.txt     # BioFVM library configuration
    ├── cmake/             # CMake helper scripts
    │   └── attach_kernels.cmake
    ├── include/           # Public API headers
    │   └── biofvm/
    │       ├── microenvironment.h       # Domain discretization
    │       ├── solver.h                 # Main solver interface
    │       ├── solver_registry.h        # Kernel registration system
    │       └── ...
    ├── src/               # Implementation files (private)
    │   ├── microenvironment.cpp
    │   ├── solver.cpp
    │   ├── solver_registry_attach.cpp   # Links all kernels
    │   └── ...
    ├── kernels/           # Hardware-specific solver backends
    │   ├── openmp_solver/
    │   │   ├── CMakeLists.txt
    │   │   ├── include/
    │   │   │   └── biofvm/kernels/openmp_solver/
    │   │   │       └── diffusion_solver.h
    │   │   └── src/
    │   │       ├── diffusion_solver.cpp
    │   │       └── register_solver.cpp  # Self-registration
    │   └── thrust_solver/
    │       ├── CMakeLists.txt
    │       ├── include/
    │       │   └── biofvm/kernels/thrust_solver/
    │       │       └── diffusion_solver.h
    │       └── src/
    │           ├── diffusion_solver.cpp
    │           └── register_solver.cpp
    ├── examples/          # Demonstration applications
    │   ├── CMakeLists.txt
    │   ├── diffuse.cpp              # Simple diffusion example
    │   └── settings.xml             # Configuration file
    └── tests/             # Unit tests
        ├── CMakeLists.txt
        ├── microenvironment_test.cpp
        └── ...
```

**Key Points:**
- `include/biofvm/` - Public API for BioFVM module
- `src/` - Private implementation details
- `kernels/` - Self-contained solver backends (OpenMP, Thrust, etc.)
- `examples/` - Sample applications demonstrating usage
- `tests/` - Unit tests using GoogleTest

### Mechanics Module

The mechanics module provides cell-cell and cell-substrate mechanical interaction models.

```
mechanics/
└── micromechanics/        # Micromechanics implementation
    ├── CMakeLists.txt     # Micromechanics library configuration
    ├── include/           # Public API headers
    │   └── micromechanics/
    │       ├── forces.h
    │       └── ...
    ├── src/               # Implementation files (private)
    │   └── forces.cpp
    └── tests/             # Unit tests
        ├── CMakeLists.txt
        └── forces_test.cpp
```

**Key Points:**
- Similar structure to reactions-diffusion
- May have `kernels/` subdirectory for alternative implementations
- Tests validate force calculations and mechanical updates

### Phenotype Module

The phenotype module integrates diffusion and mechanics into complete simulations.

```
phenotype/
└── physicore/             # PhysiCore phenotype implementation
    ├── CMakeLists.txt     # Phenotype library configuration
    └── src/               # Implementation files
        └── ...
```

**Key Points:**
- Coordinates all other modules
- Implements biological rules and cell behaviors
- Creates executable applications

## Build Output Structure

When you build PhysiCore, CMake generates output in the `build/` directory organized by preset:

```
build/
├── gcc-release/           # GCC release build
│   ├── vcpkg_installed/   # Dependencies installed by vcpkg
│   │   └── x64-linux/
│   │       ├── include/
│   │       └── lib/
│   ├── common/
│   │   ├── libphysicore.common.a
│   │   └── tests/
│   ├── reactions-diffusion/
│   │   └── biofvm/
│   │       ├── libphysicore.reactions-diffusion.biofvm.a
│   │       ├── kernels/
│   │       │   ├── openmp_solver/
│   │       │   └── thrust_solver/
│   │       ├── examples/
│   │       │   └── reactions-diffusion.biofvm.diffuse
│   │       └── tests/
│   └── ...
├── llvm-debug/            # Clang debug build with sanitizers
└── msvc-release/          # MSVC release build
```

**Key Points:**
- Each preset has its own build tree
- `vcpkg_installed/` contains dependency artifacts
- Libraries are built with `.a` (Linux/macOS) or `.lib` (Windows) extensions
- Executables in `examples/` subdirectories

## vcpkg Integration

PhysiCore uses vcpkg for dependency management with custom overlays:

```
vcpkg/                     # vcpkg submodule (git submodule)
vcpkg.json                 # Manifest: declares dependencies
vcpkg-configuration.json   # Baseline and registry configuration
ports-overlays/            # Custom or modified ports
│   └── noarr-structures/  # Example custom port
│       ├── portfile.cmake
│       └── vcpkg.json
triplets-overlays/         # Custom triplet definitions
    └── arm64-osx-gcc.cmake
```

**Key Points:**
- `vcpkg/` is a git submodule pointing to Microsoft's vcpkg repository
- `vcpkg.json` lists all dependencies (gtest, highway, tbb, etc.)
- `ports-overlays/` allows customizing dependency builds
- `triplets-overlays/` defines platform-specific build configurations

## GitHub Workflows

Continuous Integration and automation workflows:

```
.github/
├── workflows/
│   ├── cmake-multi-platform.yml     # Cross-platform build tests
│   ├── cmake-ubuntu-sanitized.yml   # Sanitizer builds (ASAN, TSAN)
│   ├── formal-checks.yml            # Linting and formatting checks
│   ├── release.yml                  # Semantic release automation
│   └── sonarqube.yml                # Code quality analysis
└── actions/
    ├── cuda-toolkit/                # Reusable CUDA setup action
    └── setup-vcpkg/                 # Reusable vcpkg setup action
```

**Key Points:**
- Multi-platform testing on Linux, macOS, Windows
- Sanitizer builds catch memory and threading issues
- Automated semantic versioning and releases

## Configuration Files

### Code Quality

- `.clang-format` - Enforces consistent code formatting
- `.clang-tidy` - Static analysis rules
- `.clangd` - Language server configuration for IDE support

### Build Configuration

- `CMakeLists.txt` - Root build configuration
- `CMakePresets.json` - Pre-defined build configurations for different compilers and platforms
- `vcpkg.json` - Dependency manifest
- `vcpkg-configuration.json` - vcpkg baseline and overlays

### Sanitizers

- `lsan.supp` - Leak Sanitizer suppressions for known false positives
- `tsan.supp` - Thread Sanitizer suppressions

## Finding Components

Use this quick reference to locate specific functionality:

| What you need | Where to find it |
|---------------|------------------|
| Core interfaces | `common/include/common/` |
| Agent data structures | `common/include/common/base_agent*.h` |
| Timestep executor | `common/include/common/timestep_executor.h` |
| BioFVM diffusion | `reactions-diffusion/biofvm/include/biofvm/` |
| OpenMP solver | `reactions-diffusion/biofvm/kernels/openmp_solver/` |
| Thrust solver | `reactions-diffusion/biofvm/kernels/thrust_solver/` |
| Mechanics forces | `mechanics/micromechanics/include/micromechanics/` |
| Example simulations | `reactions-diffusion/biofvm/examples/` |
| Unit tests | `*/tests/` in each module |
| Build presets | `CMakePresets.json` |
| Dependencies | `vcpkg.json` |

## Navigating the Code

### Reading Code
1. Start with public headers in `include/` directories
2. Check `CMakeLists.txt` to understand dependencies
3. Look at `examples/` for usage patterns
4. Refer to `tests/` for expected behavior

### Adding New Features
1. Identify the appropriate module (common, reactions-diffusion, mechanics, phenotype)
2. Add public API to `include/` directory
3. Implement in `src/` directory
4. Add tests in `tests/` directory
5. Update `CMakeLists.txt` if adding new files
6. For new solver backends, create a `kernels/` subdirectory

### Running Tests
```bash
# Build and run all tests for a preset
cmake --workflow --preset=gcc-release

# Run specific test
ctest --preset gcc-release -R microenvironment_test --output-on-failure
```

## Summary

The PhysiCore repository is organized for:
- **Clarity** - Consistent structure across all modules
- **Modularity** - Each component is self-contained
- **Discoverability** - Public APIs clearly separated from implementation
- **Testability** - Tests alongside the code they validate
- **Reproducibility** - vcpkg ensures consistent dependencies

Understanding this structure will help you efficiently navigate, modify, and extend PhysiCore for your research needs.
