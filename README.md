[![CMake Build on Ubuntu, Windows, and MacOS](https://github.com/bsc-life/PhysiCore/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/bsc-life/PhysiCore/actions/workflows/cmake-multi-platform.yml)
[![ASAN, LSAN, UBSAN, TSAN Build on Ubuntu](https://github.com/bsc-life/PhysiCore/actions/workflows/cmake-ubuntu-sanitized.yml/badge.svg)](https://github.com/bsc-life/PhysiCore/actions/workflows/cmake-ubuntu-sanitized.yml)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=bsc-life_PhysiCore&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=bsc-life_PhysiCore)

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
