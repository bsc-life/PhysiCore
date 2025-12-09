---
layout: default
title: Installation
nav_order: 2
description: "Complete installation guide for PhysiCore including dependencies, build instructions, and platform support"
---

# Installation

This guide walks you through setting up PhysiCore on your system, from installing dependencies to building your first example.

## Prerequisites

### Required Tools

- **CMake 3.24+** - Build system with presets support
- **C++20 Compiler** - See [Supported Compilers](#supported-compilers) below
- **Ninja** - Build tool (CMake presets use the Ninja generator)
- **vcpkg** - C++ package manager for dependencies

### Supported Operating Systems

PhysiCore is tested and supported on:

- **Linux** - Ubuntu 22.04+ (x86, arm64)
- **macOS** - macOS 15+ (x86, arm64)
- **Windows** - Windows 11 with Visual Studio 2019+ (x86, arm64)

### Supported Compilers

| Platform | Compiler | Minimum Version |
|----------|----------|-----------------|
| Linux | GCC | 13+ |
| Linux | Clang/LLVM | 18+ |
| macOS | GCC/AppleClang | 13+/18+ |
| Windows | MSVC/LLVM | 19.29+ (VS 2019) |

## Dependencies

PhysiCore uses vcpkg in manifest mode to manage its dependencies. The following libraries are automatically installed during the first build:

- **cccl** - NVIDIA Thrust library for TBB and CUDA backends (not on macOS)
- **gtest** - Google Test framework for unit testing
- **highway** - Vendor-agnostic SIMD vectorization library
- **noarr-structures** - Memory layouts library for efficient data access patterns
- **pugixml** - Lightweight XML parser for configuration files
- **tbb** - Intel Threading Building Blocks for parallelism (not on macOS)
- **vtk-ioxml** - VTK I/O library for serialization and visualization output

Dependencies are specified in `vcpkg.json` and automatically resolved during the CMake configure step.

## Setting Up vcpkg

PhysiCore requires vcpkg to be available via the `VCPKG_ROOT` environment variable. You have two options:

### Option A: Use the Repository's Local vcpkg (Recommended)

This ensures reproducible builds with pinned dependencies:

**Linux/macOS (bash/zsh):**
```bash
git clone --recursive https://github.com/bsc-life/PhysiCore.git
cd PhysiCore
export VCPKG_ROOT="$PWD/vcpkg"
"$VCPKG_ROOT"/bootstrap-vcpkg.sh
```

**Windows (PowerShell):**
```powershell
git clone --recursive https://github.com/bsc-life/PhysiCore.git
cd PhysiCore
$env:VCPKG_ROOT = "$PWD/vcpkg"
& "$env:VCPKG_ROOT/bootstrap-vcpkg.bat"
```

### Option B: Use an Existing vcpkg Installation

If you already have vcpkg installed:

**Linux/macOS:**
```bash
export VCPKG_ROOT="/path/to/your/vcpkg"
```

**Windows:**
```powershell
$env:VCPKG_ROOT = "C:\path\to\vcpkg"
```

> **Note:** Make sure your vcpkg installation is up to date by running `git pull` in the vcpkg directory.

## Building PhysiCore

PhysiCore uses CMake presets for a streamlined build experience. All build configurations are defined in `CMakePresets.json`.

### Quick Start: Configure, Build, and Test

**Linux/macOS (GCC):**
```bash
# Full workflow: configure + build + test
cmake --workflow --preset=gcc-release

# Or step by step:
cmake --preset=gcc-release
cmake --build --preset=gcc-release
ctest --preset gcc-release --output-on-failure
```

**Windows (MSVC):**
```powershell
cmake --workflow --preset=msvc-release
```

### Available Build Presets

The following presets are available for different compilers and build types:

#### GCC Presets
- `gcc-debug` - Debug build with ASAN, LSAN, UBSAN sanitizers
- `gcc-relwithdebinfo` - Release with debug info + TSAN
- `gcc-release` - Optimized release build

#### LLVM/Clang Presets
- `llvm-debug` - Debug build with sanitizers
- `llvm-relwithdebinfo` - Release with debug info + TSAN
- `llvm-release` - Optimized release build

#### AppleClang Presets (macOS)
- `appleclang-debug` - Debug build with sanitizers
- `appleclang-relwithdebinfo` - Release with debug info + TSAN
- `appleclang-release` - Optimized release build

#### MSVC Presets (Windows)
- `msvc-debug` - Debug build
- `msvc-release` - Optimized release build

### Build Configuration Details

On the first configure, vcpkg will:
1. Download and build all dependencies
2. Install them under `build/<preset>/vcpkg_installed/`
3. Generate build files for the project

Subsequent builds reuse the cached dependencies unless `vcpkg.json` changes.

## Installing PhysiCore

To install PhysiCore libraries and headers to a system location:

```bash
cmake --preset=gcc-release
cmake --build --preset=gcc-release
cmake --install build/gcc-release --prefix /usr/local
```

This installs:
- Libraries to `<prefix>/lib/`
- Headers to `<prefix>/include/`
- CMake config files to `<prefix>/lib/cmake/physicore/`

## Running an Example

PhysiCore includes example applications to demonstrate its capabilities. Here's how to build and run the diffusion example:

### Build the Diffusion Example

```bash
cmake --preset=gcc-release
cmake --build --preset=gcc-release --target reactions-diffusion.biofvm.diffuse
```

### Run the Example

```bash
./build/gcc-release/reactions-diffusion/biofvm/examples/reactions-diffusion.biofvm.diffuse
```

The simulation runs for 30 minutes of simulated time and generates VTK output files in the `output/` directory for visualization with tools like ParaView.

### Visualizing Results

The output VTK files can be opened in [ParaView](https://www.paraview.org/) or other VTK-compatible visualization tools:

```bash
paraview output/*.vtu
```

## Troubleshooting

### vcpkg Not Found

If you see `CMAKE_TOOLCHAIN_FILE not found` errors:
- Ensure `VCPKG_ROOT` is set correctly
- Verify the path exists: `ls $VCPKG_ROOT` (Linux/macOS) or `dir $env:VCPKG_ROOT` (Windows)
- Bootstrap vcpkg if you haven't already

### Compiler Not Found

If CMake can't find your compiler:
- Install the required compiler (see [Supported Compilers](#supported-compilers))
- Ensure it's in your system PATH
- On Linux, install with: `sudo apt install gcc g++` or `sudo apt install clang`

### CUDA Support

CUDA support is optional and automatically detected when available. The standard presets (e.g., `gcc-release`, `llvm-release`) will use CUDA if the NVIDIA CUDA Toolkit is installed.

To enable CUDA support:
1. Install [NVIDIA CUDA Toolkit](https://developer.nvidia.com/cuda-downloads)
2. Use any standard preset (e.g., `gcc-release`, `llvm-release`)
3. CMake will automatically detect and enable CUDA support

## Next Steps

- Explore the [Architecture](Architecture.md) to understand how PhysiCore is organized
- Browse the [Repository Structure](Repository-Structure.md) to navigate the codebase
