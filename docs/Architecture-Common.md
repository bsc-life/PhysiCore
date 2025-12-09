---
layout: default
title: Common Module
parent: Architecture
nav_order: 1
description: "Core abstractions, timestep executors, agent containers, and foundational types in PhysiCore"
---

# Common Module
{: .no_toc }

The `common` module provides the foundational abstractions that all other PhysiCore modules build upon. It defines core interfaces, data structures, and types that enable modular simulation design.

## Table of Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Overview

The Common module establishes the contracts and base types that enable PhysiCore's modular architecture. Every other module depends on Common, making it the foundation of the framework.

**Namespace:** `physicore::common`

**Location:** `common/include/common/`

**Key Responsibilities:**
- Define the simulation loop contract via `timestep_executor`
- Provide agent data structures using structure-of-arrays (SoA) pattern
- Manage agent collections through containers
- Establish core types and concepts for type safety

## The `timestep_executor` Interface

The `timestep_executor` is the central abstraction that defines how simulation components participate in the main simulation loop.

### Interface Definition

```cpp
namespace physicore::common {

class timestep_executor {
public:
    virtual ~timestep_executor() = default;
    
    // Execute a single timestep of the simulation
    virtual void run_single_timestep() = 0;
    
    // Serialize current simulation state
    virtual void serialize_state(real_t current_time) = 0;
};

} // namespace physicore::common
```

### Design Rationale

All simulation components (diffusion solvers, mechanics engines, phenotype models) implement this interface, enabling:

1. **Uniform simulation loop** - The main loop doesn't need to know implementation details
2. **Composability** - Multiple executors can be orchestrated together
3. **Testability** - Each component can be tested in isolation
4. **Extensibility** - New modules automatically integrate with existing infrastructure

### Usage Example

```cpp
#include <common/timestep_executor.h>

class DiffusionSolver : public physicore::common::timestep_executor {
public:
    void run_single_timestep() override {
        // Solve reaction-diffusion PDE for one timestep
        compute_diffusion();
        apply_reactions();
        update_concentrations();
    }
    
    void serialize_state(real_t current_time) override {
        // Write substrate concentrations to VTK files
        vtk_writer.write(current_time, substrate_data);
    }
};
```

## Agent Data Structures

PhysiCore uses a **structure-of-arrays (SoA)** pattern for agent data to enable efficient vectorization and cache-friendly memory access.

### Structure-of-Arrays (SoA) Pattern

Instead of storing agent data as an array of structures (AoS):

```cpp
// Array-of-Structures (AoS) - NOT used in PhysiCore
struct Agent {
    double x, y, z;      // Position
    double vx, vy, vz;   // Velocity
    double radius;
    int type;
};
std::vector<Agent> agents;  // Poor cache locality
```

PhysiCore uses structure-of-arrays (SoA):

```cpp
// Structure-of-Arrays (SoA) - Used in PhysiCore
struct AgentData {
    std::vector<double> x, y, z;        // Positions
    std::vector<double> vx, vy, vz;     // Velocities
    std::vector<double> radius;
    std::vector<int> type;
};
```

**Benefits:**
- **Vectorization** - SIMD operations process contiguous data efficiently
- **Cache efficiency** - Related data is stored together
- **Memory bandwidth** - Fewer cache misses during computation
- **Parallelization** - Easy to partition across threads/GPUs

### The `base_agent_data` Class

The `base_agent_data` class provides the foundational SoA storage:

```cpp
namespace physicore::common {

class base_agent_data {
public:
    // Positions
    std::vector<real_t> positions_x;
    std::vector<real_t> positions_y;
    std::vector<real_t> positions_z;
    
    // Velocities
    std::vector<real_t> velocities_x;
    std::vector<real_t> velocities_y;
    std::vector<real_t> velocities_z;
    
    // Add a new agent
    void push_back(real_t x, real_t y, real_t z,
                   real_t vx, real_t vy, real_t vz);
    
    // Get number of agents
    std::size_t size() const;
    
    // Remove agent at index
    void erase(std::size_t index);
};

} // namespace physicore::common
```

### The `base_agent` Proxy

To provide object-like access to individual agents while maintaining SoA storage, PhysiCore uses a **proxy pattern** via the `base_agent` class:

```cpp
namespace physicore::common {

class base_agent {
public:
    // Constructor takes reference to data and index
    base_agent(base_agent_data& data, std::size_t index);
    
    // Access position
    real_t position_x() const;
    real_t position_y() const;
    real_t position_z() const;
    
    // Modify position
    void set_position(real_t x, real_t y, real_t z);
    
    // Access velocity
    real_t velocity_x() const;
    void set_velocity_x(real_t vx);
    
    // ... similar for other properties
};

} // namespace physicore::common
```

**Usage:**

```cpp
base_agent_data agents;
agents.push_back(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);

// Access via proxy
base_agent agent(agents, 0);
agent.set_position(1.0, 2.0, 3.0);

real_t x = agent.position_x();  // Returns 1.0
```

## Agent Containers

The `base_agent_container` provides high-level management of agent collections:

```cpp
namespace physicore::common {

class base_agent_container {
public:
    // Add a new agent
    base_agent create_agent(real_t x, real_t y, real_t z);
    
    // Access agent by index
    base_agent operator[](std::size_t index);
    
    // Iterate over all agents
    auto begin();
    auto end();
    
    // Get number of agents
    std::size_t size() const;
    
    // Remove agent
    void remove(std::size_t index);
    
private:
    base_agent_data data_;
};

} // namespace physicore::common
```

**Usage Example:**

```cpp
#include <common/base_agent_container.h>

base_agent_container container;

// Create agents
auto agent1 = container.create_agent(0.0, 0.0, 0.0);
auto agent2 = container.create_agent(1.0, 1.0, 1.0);

// Iterate and modify
for (auto agent : container) {
    real_t x = agent.position_x();
    agent.set_velocity_x(x * 0.1);
}
```

## Core Types

The Common module defines fundamental types used throughout PhysiCore:

```cpp
namespace physicore::common {

// Floating-point precision (configurable)
using real_t = double;

// Index types
using index_t = std::size_t;

// Time types
using time_t = real_t;

} // namespace physicore::common
```

## C++20 Concepts

PhysiCore uses C++20 concepts to enforce type constraints and provide better error messages:

```cpp
namespace physicore::common {

// Concept: Type must be a timestep executor
template<typename T>
concept TimestepExecutor = requires(T t, real_t time) {
    { t.run_single_timestep() } -> std::same_as<void>;
    { t.serialize_state(time) } -> std::same_as<void>;
};

// Concept: Type must be an agent container
template<typename T>
concept AgentContainer = requires(T t, std::size_t index) {
    { t.size() } -> std::convertible_to<std::size_t>;
    { t[index] };
    { t.begin() };
    { t.end() };
};

} // namespace physicore::common
```

## File Organization

The Common module follows PhysiCore's public/private API separation:

```
common/
├── include/
│   └── common/              # Public API headers
│       ├── timestep_executor.h
│       ├── base_agent.h
│       ├── base_agent_data.h
│       ├── base_agent_container.h
│       ├── base_agent_interface.h
│       ├── generic_agent_container.h
│       ├── generic_agent_solver.h
│       ├── concepts.h
│       └── types.h
├── src/                     # Private implementation (if needed)
└── tests/
    ├── test_base_agent.cpp
    ├── test_base_agent_data.cpp
    └── test_base_agent_container.cpp
```

## Public API Stability

Headers in `common/include/common/` are exported via CMake `FILE_SET HEADERS` and constitute the **stable public API**. These interfaces are maintained across minor versions with semantic versioning guarantees.

## Dependencies

The Common module has minimal external dependencies:
- C++20 standard library
- No external libraries required

This makes it the lightest-weight component and suitable as a foundation for all other modules.

## Next Steps

Learn how the Common module interfaces are used in specific physics domains:

- **[Reactions-Diffusion Module](Architecture-Diffusion.md)** - Uses `timestep_executor` for diffusion solvers
- **[Mechanics Module](Architecture-Mechanics.md)** - Uses agent containers for force calculations
- **[Phenotype Module](Architecture-Phenotype.md)** - Orchestrates multiple executors

---

**See also:**
- [Architecture Overview](Architecture.md)
- [Build System Integration](Architecture.md#build-system-integration)
