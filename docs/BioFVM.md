---
layout: default
title: BioFVM Module
nav_order: 5
description: "Complete guide to the BioFVM reaction-diffusion module including mathematical model and API reference"
has_children: false
---

# BioFVM Module
{: .no_toc }

BioFVM is PhysiCore's reaction-diffusion module for simulating substrate transport and reaction kinetics in 3D domains. It implements a finite volume method for solving partial differential equations on Cartesian meshes.

## Table of Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Mathematical Model

### Governing Equations

BioFVM solves the reaction-diffusion-decay equation for multiple substrates in 3D space:

$$
\frac{\partial \rho_i}{\partial t} = D_i \nabla^2 \rho_i - \lambda_i \rho_i + S_i(\mathbf{x}, t)
$$

Where:
- $$\rho_i(\mathbf{x}, t)$$ is the density (concentration) of substrate $$i$$ at position $$\mathbf{x}$$ and time $$t$$
- $$D_i$$ is the diffusion coefficient for substrate $$i$$
- $$\lambda_i$$ is the decay rate for substrate $$i$$
- $$S_i(\mathbf{x}, t)$$ represents source and sink terms (e.g., cell secretion and uptake)
- $$\nabla^2$$ is the Laplacian operator

### Discretization

The domain is discretized into a Cartesian mesh of rectangular voxels. The finite volume method is applied with:

- **Spatial discretization**: 7-point stencil for the Laplacian in 3D
- **Time integration**: First-order explicit Euler scheme
- **Boundary conditions**: Dirichlet conditions on domain boundaries and interior voxels

### Cell-Substrate Interactions

Cells contribute to the source/sink term $$S_i$$ through:

1. **Secretion**: Cells release substrates into the microenvironment
2. **Uptake**: Cells consume substrates from their local voxel
3. **Internalization**: Optional tracking of substrate amounts inside cells

The uptake rate can follow Michaelis-Menten kinetics or linear models depending on configuration.

---

## Core API

### Microenvironment Class

The `microenvironment` class is the main entry point for BioFVM simulations.

#### Construction

```cpp
#include <biofvm/microenvironment.h>

using namespace physicore::biofvm;

// Create a Cartesian mesh
cartesian_mesh mesh(
    3,                                    // 3D domain
    {-500, -500, -500},                   // Bounding box min
    {500, 500, 500},                      // Bounding box max
    {20, 20, 20}                          // Voxel dimensions
);

// Create microenvironment with 2 substrates
auto m = std::make_unique<microenvironment>(
    mesh,
    2,                                    // Number of substrates
    0.01                                  // Timestep (minutes)
);
```

#### From Configuration File

```cpp
// Load from XML configuration
auto m = microenvironment::create_from_config("settings.xml");
```

#### Key Methods

##### Run Simulation

```cpp
// Execute a single timestep
m->run_single_timestep();

// The simulation time is automatically advanced
std::cout << "Time: " << m->simulation_time << std::endl;
```

##### Query Substrate Densities

```cpp
// Get substrate density at a specific voxel
real_t density = m->get_substrate_density(
    substrate_idx,  // Substrate index
    x, y, z         // Voxel coordinates
);
```

##### Dirichlet Boundary Conditions

```cpp
// Set Dirichlet condition on minimum X boundary
m->update_dirichlet_boundary_min(
    'x',                // Dimension ('x', 'y', or 'z')
    substrate_idx,      // Substrate index
    21.0,               // Concentration value
    true                // Enable condition
);

// Set Dirichlet condition on maximum Y boundary
m->update_dirichlet_boundary_max(
    'y',
    substrate_idx,
    0.0,
    true
);

// Set interior Dirichlet voxel (e.g., blood vessel)
m->update_dirichlet_interior_voxel(
    {50, 50, 50},       // Voxel position
    substrate_idx,
    38.0,               // Oxygen concentration
    true
);

// Apply all Dirichlet condition updates
m->update_dirichlet_conditions();
```

##### Serialization

```cpp
// Save current state to VTK files
m->serialize_state(current_time);
```

#### Configuration Parameters

The `microenvironment` object exposes configuration through public members:

```cpp
// Environment metadata
m->name = "tumor_microenvironment";
m->time_units = "min";
m->space_units = "micron";

// Substrate properties
m->substrates_names[0] = "oxygen";
m->substrates_names[1] = "glucose";
m->substrates_units[0] = "mmHg";
m->substrates_units[1] = "mM";

// Diffusion and decay parameters
m->diffusion_coefficients[0] = 100000.0;  // micron²/min
m->decay_rates[0] = 0.1;                  // 1/min

m->diffusion_coefficients[1] = 50000.0;
m->decay_rates[1] = 0.01;

// Initial conditions
m->initial_conditions[0] = 38.0;  // mmHg
m->initial_conditions[1] = 5.0;   // mM
```

---

### Cartesian Mesh

The `cartesian_mesh` struct defines the spatial discretization.

#### Construction

```cpp
cartesian_mesh mesh(
    dims,                 // 1, 2, or 3 dimensions
    bounding_box_mins,    // {x_min, y_min, z_min}
    bounding_box_maxs,    // {x_max, y_max, z_max}
    voxel_shape          // {dx, dy, dz}
);
```

#### Key Methods

```cpp
// Get total number of voxels
std::size_t count = mesh.voxel_count();

// Get volume of a single voxel
index_t vol = mesh.voxel_volume();

// Convert 3D indices to linear index
std::size_t idx = mesh.linearize(x, y, z);

// Find voxel containing a position
std::array<real_t, 3> position = {100.0, -50.0, 25.0};
auto [vx, vy, vz] = mesh.voxel_position(position);

// Get center coordinates of a voxel
auto center = mesh.voxel_center({vx, vy, vz});
```

---

### Solver Interface

BioFVM supports multiple solver backends through the `solver` interface.

#### Solver Selection

```cpp
#include <biofvm/solver_registry.h>

// Available solvers: "openmp", "thrust"
auto solver = solver_registry::instance().get("openmp");
m->solver = std::move(solver);
```

#### Solver Methods

```cpp
class solver {
public:
    // Initialize substrate densities
    virtual void initialize(microenvironment& m) = 0;
    
    // Solve diffusion-decay for N iterations
    virtual void solve(microenvironment& m, index_t iterations) = 0;
    
    // Access substrate density
    virtual real_t get_substrate_density(
        index_t s, index_t x, index_t y, index_t z
    ) const = 0;
    
    // Update Dirichlet conditions
    virtual void reinitialize_dirichlet(microenvironment& m) = 0;
    
    // GPU solvers: transfer data
    virtual void transfer_to_device(microenvironment& m);
    virtual void transfer_to_host(microenvironment& m);
};
```

---

### Agent Container

Agents (cells) are managed through the agent container interface.

#### Agent Interface

```cpp
#include <biofvm/agent.h>

class agent_interface {
public:
    // Position and volume
    virtual std::span<const real_t> position() const = 0;
    virtual real_t volume() const = 0;
    
    // Substrate interactions
    virtual real_t secretion_rate(index_t substrate_index) const = 0;
    virtual real_t uptake_rate(index_t substrate_index) const = 0;
    virtual real_t saturation_density(index_t substrate_index) const = 0;
    
    // Internalized substrates (optional)
    virtual real_t& internalized_substrates(index_t substrate_index) = 0;
};
```

#### Container Access

```cpp
// Access agents container
container_ptr agents = m->agents;

// Iterate over agents
for (auto& agent : *agents) {
    auto pos = agent.position();
    auto vol = agent.volume();
    
    // Get uptake rate for oxygen
    auto uptake = agent.uptake_rate(0);
}
```

---

## Solver Backends

### OpenMP Solver

CPU-parallel solver using OpenMP for shared-memory parallelism.

**Features:**
- Multi-threaded diffusion computation
- Efficient for CPU-only systems
- Automatic thread scaling

**Usage:**
```cpp
auto solver = solver_registry::instance().get("openmp");
```

### Thrust Solver

GPU/CPU solver using NVIDIA Thrust library with TBB backend.

**Features:**
- Supports both CPU (via TBB) and GPU (via CUDA) execution
- Vectorized operations
- Efficient for large-scale simulations

**Usage:**
```cpp
auto solver = solver_registry::instance().get("thrust");

// For GPU execution, transfer data before solving
solver->transfer_to_device(*m);
solver->solve(*m, iterations);
solver->transfer_to_host(*m);
```

---

## Complete Example

Here's a complete example demonstrating BioFVM usage:

```cpp
#include <biofvm/microenvironment.h>
#include <biofvm/solver_registry.h>
#include <iostream>

using namespace physicore::biofvm;

int main() {
    // Create 3D mesh: 1mm³ domain with 20μm voxels
    cartesian_mesh mesh(
        3,
        {-500, -500, -500},
        {500, 500, 500},
        {20, 20, 20}
    );
    
    // Create microenvironment with oxygen
    auto m = std::make_unique<microenvironment>(mesh, 1, 0.01);
    
    // Configure oxygen
    m->substrates_names[0] = "oxygen";
    m->substrates_units[0] = "mmHg";
    m->diffusion_coefficients[0] = 100000.0;  // μm²/min
    m->decay_rates[0] = 0.1;                  // 1/min
    m->initial_conditions[0] = 38.0;          // mmHg
    
    // Set boundary: normoxic on all sides
    for (char dim : {'x', 'y', 'z'}) {
        m->update_dirichlet_boundary_min(dim, 0, 38.0, true);
        m->update_dirichlet_boundary_max(dim, 0, 38.0, true);
    }
    m->update_dirichlet_conditions();
    
    // Select OpenMP solver
    m->solver = solver_registry::instance().get("openmp");
    
    // Initialize
    m->solver->initialize(*m);
    
    // Run simulation for 30 minutes
    real_t end_time = 30.0;
    while (m->simulation_time < end_time) {
        m->run_single_timestep();
        
        if (static_cast<int>(m->simulation_time) % 10 == 0) {
            std::cout << "Time: " << m->simulation_time 
                      << " min" << std::endl;
            m->serialize_state(m->simulation_time);
        }
    }
    
    std::cout << "Simulation complete!" << std::endl;
    return 0;
}
```

---

## Advanced Topics

### Custom Bulk Sources/Sinks

Implement custom spatially-varying source terms:

```cpp
#include <biofvm/bulk_functor.h>

class CustomBulk : public bulk_functor {
public:
    real_t evaluate(
        index_t substrate_idx,
        index_t voxel_idx,
        const microenvironment& m
    ) const override {
        // Example: Gaussian source at origin
        // First, convert linear index to 3D coordinates
        index_t z = voxel_idx / (m.mesh.grid_shape[0] * m.mesh.grid_shape[1]);
        index_t rem = voxel_idx % (m.mesh.grid_shape[0] * m.mesh.grid_shape[1]);
        index_t y = rem / m.mesh.grid_shape[0];
        index_t x = rem % m.mesh.grid_shape[0];
        
        // Get voxel center position
        auto center = m.mesh.voxel_center({x, y, z});
        real_t r2 = center[0]*center[0] + 
                    center[1]*center[1] + 
                    center[2]*center[2];
        return 1000.0 * std::exp(-r2 / 10000.0);
    }
};

m->bulk_fnc = std::make_unique<CustomBulk>();
```

### Performance Tuning

**Timestep Selection:**
- Choose timestep based on CFL condition: $$\Delta t \leq \frac{\Delta x^2}{6D}$$
- Smaller timesteps increase accuracy but require more iterations

**Voxel Resolution:**
- Finer voxels capture sharper gradients but increase computational cost
- Balance resolution with available memory and compute resources

**Solver Choice:**
- OpenMP: Best for moderate-size CPU-only simulations
- Thrust/TBB: Better for very large CPU simulations with vectorization
- Thrust/CUDA: Essential for GPU-accelerated large-scale simulations

---

## Further Reading

- [Architecture Guide](Architecture.md) - Understanding PhysiCore's module system
- [Installation](Installation.md) - Building with BioFVM support
- [Repository Structure](Repository-Structure.md) - Navigating BioFVM source code
- [Original BioFVM Paper](http://dx.doi.org/10.1093/bioinformatics/btv730) - Mathematical foundations
