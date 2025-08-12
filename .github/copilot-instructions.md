# Copilot Instructions for PhysiCore

## Project Overview
- **PhysiCore** is a modular, high-performance agent-based simulation framework, re-architected from PhysiCell.
- Major components:
  - `common/`: Core agent data structures and utilities (e.g., `base_agent_data.h`, `types.h`).
  - `mechanics/physicell/`: Mechanics engines (force models, integration methods).
  - `reactions-diffusion/biofvm/`: Diffusion solvers (FEM, FVM, FDM, etc.).
  - `src/`: Main entry point and top-level simulation logic.
- Each submodule (mechanics, diffusion) is designed to be swappable and independently testable.

## Build & Test Workflow
- **Build:**
  - Use CMake (see `CMakeLists.txt` at root and in submodules).
  - Typical build: `cmake -S . -B build && cmake --build build`
- **Test:**
  - Tests are in `common/src/tests/`, `mechanics/physicell/tests/`, `reactions-diffusion/biofvm/tests/`.
  - Run all tests: `ctest --test-dir build`

## Coding Patterns & Conventions
- **Namespaces:** All core code is under the `physicore` namespace.
- **Data Layout:**
  - Agent data is stored in SoA (Structure of Arrays) style for performance (see `base_agent_data`).
  - Use `index_t` and `real_t` typedefs from `types.h` for indices and floating-point values.
- **Extensibility:**
  - Add new mechanics or diffusion models by creating new subfolders under `mechanics/` or `reactions-diffusion/` and updating CMake.
  - Keep interfaces minimal and decoupled; avoid cross-module dependencies.
- **Testing:**
  - Use GoogleTest (vendored in `build/_deps/googletest-src/`).
  - Place new tests in the relevant `tests/` subfolder.

## Integration & External Dependencies
- **GoogleTest** is included as a submodule and built automatically.
- No external runtime dependencies beyond standard C++ and CMake.

## Examples
- See `common/include/base_agent_data.h` for agent data layout and API.
- See `mechanics/physicell/` and `reactions-diffusion/biofvm/` for pluggable module structure.

## Tips for AI Agents
- When adding new modules, mirror the structure and CMake patterns of existing ones.
- Prefer explicit, minimal interfaces between modules.
- Always update or add tests for new features.
- Reference the root `README.md` for project vision and architectural goals.
