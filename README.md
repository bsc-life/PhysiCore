[![CMake Build on Ubuntu, Windows, and MacOS](https://github.com/bsc-life/PhysiCore/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/bsc-life/PhysiCore/actions/workflows/cmake-multi-platform.yml)
[![ASAN, LSAN, UBSAN, TSAN Build on Ubuntu](https://github.com/bsc-life/PhysiCore/actions/workflows/cmake-ubuntu-sanitized.yml/badge.svg)](https://github.com/bsc-life/PhysiCore/actions/workflows/cmake-ubuntu-sanitized.yml)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=bsc-life_PhysiCore&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=bsc-life_PhysiCore)
[![Release](https://img.shields.io/github/v/release/bsc-life/PhysiCore)](https://github.com/bsc-life/PhysiCore/releases/latest)
[![semantic-release: angular](https://img.shields.io/badge/semantic--release-angular-e10079?logo=semantic-release)](https://github.com/semantic-release/semantic-release)

# PhysiCore

**PhysiCore** is a modern C++20 re-architecture of the [PhysiCell](http://physicell.org) agent-based multicellular simulation framework. Built for flexibility, performance, and clarity, PhysiCore enables researchers to experiment with different numerical methods and physical models without overhauling the entire simulation engine.

## Key Features

- **Modular architecture** – Clean separation between diffusion, mechanics, and phenotype modules
- **Pluggable solvers** – Runtime-selectable backends (OpenMP, Thrust/TBB, CUDA)
- **Modern C++20** – Leveraging the latest language features for safety and performance
- **HPC-ready** – Designed for vectorization, multi-threading, and GPU acceleration
- **Cross-platform** – Tested on Linux, macOS, and Windows

## Quick Start

```bash
# Clone with submodules
git clone --recursive https://github.com/bsc-life/PhysiCore.git
cd PhysiCore

# Set up vcpkg
export VCPKG_ROOT="$PWD/vcpkg"
"$VCPKG_ROOT"/bootstrap-vcpkg.sh

# Build and test
cmake --workflow --preset=gcc-release

# Run an example
./build/gcc-release/reactions-diffusion/biofvm/examples/reactions-diffusion.biofvm.diffuse
```

## Documentation

Comprehensive documentation is available in the [`docs/`](docs/) directory:

- **[Installation Guide](docs/Installation.md)** - Dependencies, build instructions, and platform support
- **[Architecture](docs/Architecture.md)** - Modular design, implementations, and solver backends
- **[Repository Structure](docs/Repository-Structure.md)** - Directory layout and code organization

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines and [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) for community standards.

## License

PhysiCore is released under the BSD 3-Clause License. See [LICENSE](LICENSE) for details.
