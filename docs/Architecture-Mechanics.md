---
layout: default
title: Mechanics Module
parent: Architecture
nav_order: 3
description: "Cell-cell and cell-substrate mechanical interactions using force-based models in PhysiCore"
---

# Mechanics Module
{: .no_toc }

The mechanics module handles cell-cell and cell-substrate mechanical interactions in multicellular simulations. It computes forces between cells and updates cell positions based on spring-based models, repulsion, and adhesion.

## Table of Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Overview

Mechanical interactions operate on timescales of **minutes to hours**, slower than substrate diffusion but faster than cell division. The mechanics module provides:

- **Cell-cell repulsion** - Prevents cell overlap
- **Cell-cell adhesion** - Holds tissues together
- **Spring forces** - Elastic cell deformation
- **Boundary forces** - Domain confinement
- **Position updates** - Integration of equations of motion

**Namespace:** `physicore::mechanics`

**Location:** `mechanics/`

## Physical Principles

### Force Balance

Each cell experiences multiple forces that determine its motion:

$$
m_i \frac{d^2\mathbf{x}_i}{dt^2} = \mathbf{F}_i^{\text{repulsion}} + \mathbf{F}_i^{\text{adhesion}} + \mathbf{F}_i^{\text{spring}} + \mathbf{F}_i^{\text{boundary}} + \mathbf{F}_i^{\text{drag}}
$$

Where:
- $\mathbf{x}_i$ - Position of cell $i$
- $m_i$ - Cell mass
- $\mathbf{F}_i$ - Various force contributions

### Overdamped Dynamics

At cellular scales, inertia is negligible compared to viscous drag. The system operates in the **overdamped regime**:

$$
\gamma \frac{d\mathbf{x}_i}{dt} = \sum \mathbf{F}_i
$$

Where $\gamma$ is the drag coefficient. This simplifies to:

$$
\frac{d\mathbf{x}_i}{dt} = \frac{1}{\gamma} \sum \mathbf{F}_i
$$

**Physical interpretation:**
- Cells reach terminal velocity instantly
- No momentum or acceleration
- Position updates are first-order in forces

## Force Models

### Cell-Cell Repulsion

When cells overlap, they experience repulsive forces to prevent interpenetration:

$$
\mathbf{F}_{ij}^{\text{repulsion}} = 
\begin{cases}
k_{\text{rep}} (d_{ij} - r_i - r_j) \frac{\mathbf{x}_j - \mathbf{x}_i}{d_{ij}} & \text{if } d_{ij} < r_i + r_j \\
\mathbf{0} & \text{otherwise}
\end{cases}
$$

Where:
- $d_{ij} = |\mathbf{x}_j - \mathbf{x}_i|$ - Distance between cell centers
- $r_i, r_j$ - Cell radii
- $k_{\text{rep}}$ - Repulsion coefficient (stiffness)
- $d_{ij} < r_i + r_j$ - Cells are overlapping

**Properties:**
- **Magnitude** proportional to overlap distance
- **Direction** pushes cells apart
- **Range** only active when overlapping
- **Stiffness** determines how "hard" cells are

### Cell-Cell Adhesion

Cells within adhesion range attract each other:

$$
\mathbf{F}_{ij}^{\text{adhesion}} = 
\begin{cases}
-k_{\text{adh}} (d_{ij} - r_i - r_j) \frac{\mathbf{x}_j - \mathbf{x}_i}{d_{ij}} & \text{if } r_i + r_j < d_{ij} < r_{\text{max}} \\
\mathbf{0} & \text{otherwise}
\end{cases}
$$

Where:
- $k_{\text{adh}}$ - Adhesion coefficient
- $r_{\text{max}}$ - Maximum adhesion range

**Properties:**
- **Magnitude** proportional to separation distance
- **Direction** pulls cells together
- **Range** active between touching and $r_{\text{max}}$
- **Selectivity** can be cell-type specific

### Spring Forces

Cells connected by explicit springs (e.g., mother-daughter pairs):

$$
\mathbf{F}_{ij}^{\text{spring}} = -k_{\text{spring}} (d_{ij} - \ell_0) \frac{\mathbf{x}_j - \mathbf{x}_i}{d_{ij}}
$$

Where:
- $k_{\text{spring}}$ - Spring constant
- $\ell_0$ - Rest length

**Hooke's Law:**
- Tension when $d_{ij} > \ell_0$
- Compression when $d_{ij} < \ell_0$

### Boundary Forces

Confine cells within simulation domain:

$$
\mathbf{F}_i^{\text{boundary}} = -k_{\text{wall}} \max(0, r_i - x_i + x_{\text{min}}) \hat{\mathbf{x}} + \ldots
$$

Similar terms for all six domain boundaries.

## Module Interface Contract

All implementations of the mechanics module must implement the `timestep_executor` interface:

```cpp
#include <common/timestep_executor.h>

namespace physicore::mechanics {

class mechanics_solver : public physicore::common::timestep_executor {
public:
    // Run one timestep of mechanics
    void run_single_timestep() override;
    
    // Serialize cell positions
    void serialize_state(real_t current_time) override;
    
    // Module-specific interface
    virtual void compute_forces() = 0;
    virtual void update_positions(real_t dt) = 0;
    virtual void set_repulsion_coefficient(real_t k_rep) = 0;
    virtual void set_adhesion_coefficient(real_t k_adh) = 0;
};

} // namespace physicore::mechanics
```

## Implementations

### Micromechanics Implementation

**Micromechanics** is PhysiCore's spring-based mechanics engine.

**Namespace:** `physicore::mechanics::micromechanics`

**Location:** `mechanics/micromechanics/`

**Key Features:**
- Spring-based cell-cell interactions
- Cell-cell repulsion when overlapping
- Cell-cell adhesion forces
- Boundary force models
- Efficient neighbor search

**Status:** Under active development

### Force Calculation Algorithm

```cpp
void micromechanics::compute_forces() {
    // 1. Clear forces from previous timestep
    forces.assign(num_cells, {0, 0, 0});
    
    // 2. Find neighboring cells (spatial acceleration)
    auto neighbors = find_neighbors();
    
    // 3. Compute pairwise forces
    for (auto [i, j] : neighbors) {
        auto F_rep = compute_repulsion(i, j);
        auto F_adh = compute_adhesion(i, j);
        
        forces[i] += F_rep + F_adh;
        forces[j] -= F_rep + F_adh;  // Newton's 3rd law
    }
    
    // 4. Add boundary forces
    for (size_t i = 0; i < num_cells; ++i) {
        forces[i] += compute_boundary_force(i);
    }
}
```

### Position Update Algorithm

Using forward Euler integration:

```cpp
void micromechanics::update_positions(real_t dt) {
    for (size_t i = 0; i < num_cells; ++i) {
        // Velocity from force balance (overdamped)
        auto velocity = forces[i] / drag_coefficient;
        
        // Update position
        positions[i] += velocity * dt;
    }
}
```

## Neighbor Search

Finding nearby cells is the computational bottleneck. Efficient implementations use:

### Spatial Hashing
- Partition space into bins
- Hash cell positions to bins
- Check only nearby bins
- **Complexity:** $\mathcal{O}(N)$ average case

### Verlet Lists
- Cache neighbor lists
- Rebuild only when cells move significantly
- Reduces redundant distance calculations
- **Speedup:** 2-5× for typical systems

### Example: Spatial Hash

```cpp
class SpatialHash {
    std::unordered_map<int3, std::vector<size_t>> bins;
    real_t bin_size;
    
public:
    void insert(size_t cell_id, vec3 position) {
        int3 bin = floor(position / bin_size);
        bins[bin].push_back(cell_id);
    }
    
    std::vector<size_t> find_neighbors(vec3 position, real_t radius) {
        std::vector<size_t> neighbors;
        int3 center = floor(position / bin_size);
        
        // Check 3³ = 27 neighboring bins
        for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz) {
            int3 bin = center + int3{dx, dy, dz};
            if (bins.count(bin)) {
                neighbors.insert(neighbors.end(), 
                    bins[bin].begin(), bins[bin].end());
            }
        }
        return neighbors;
    }
};
```

## Time Integration

### Timestep Selection

The mechanics timestep is constrained by:

$$
\Delta t \leq \frac{\gamma}{k_{\text{max}}}
$$

Where $k_{\text{max}}$ is the largest force coefficient.

**Typical values:**
- $\gamma \sim 1$ (normalized units)
- $k_{\text{rep}} \sim 10$
- $\Delta t \leq 0.1$ (dimensionless time units)

### Integration Schemes

**Forward Euler** (used in PhysiCore):
$$
\mathbf{x}_{i}^{n+1} = \mathbf{x}_i^n + \frac{\Delta t}{\gamma} \mathbf{F}_i^n
$$

**Advantages:**
- Simple and fast
- Sufficient for overdamped systems
- Low memory overhead

**Alternatives** (future):
- Runge-Kutta 4th order - higher accuracy
- Velocity Verlet - better energy conservation

## Cell Growth and Division

Mechanics handles cell size changes:

### Growth
Cells increase radius over time:
$$
\frac{dr_i}{dt} = r_{\text{growth}}
$$

This gradually increases repulsion forces with neighbors.

### Division
When $r_i > r_{\text{divide}}$:
1. Create daughter cell at offset position
2. Reduce both radii: $r_{\text{daughter}} = r_{\text{mother}} = r_i / 2^{1/3}$
3. Optionally add spring between mother and daughter

## Coupling with Other Modules

### Reactions-Diffusion Module
- **Input:** Substrate concentrations affect cell behaviors
- **Output:** Cell positions determine source/sink locations
- **Frequency:** Mechanics updates slower than diffusion

### Phenotype Module
- **Input:** Phenotype determines force coefficients
- **Output:** Mechanical stress triggers phenotype transitions
- **Coupling:** Bidirectional feedback

See [Phenotype Module](Architecture-Phenotype.md#module-communication) for orchestration.

## Performance Optimization

### Parallelization Strategies

**Thread-level parallelism:**
```cpp
#pragma omp parallel for
for (size_t i = 0; i < num_cells; ++i) {
    forces[i] = compute_total_force(i);
}
```

**GPU acceleration:**
- Cells processed in parallel on GPU
- Spatial hashing on device memory
- Force reductions using atomics

### Memory Access Patterns

Use structure-of-arrays (SoA) from Common module:
```cpp
// Efficient: sequential access
for (size_t i = 0; i < num_cells; ++i) {
    force_x[i] = ...;  // Cache-friendly
}

// Inefficient: scattered access
for (auto& cell : cells) {
    cell.force.x = ...;  // Cache misses
}
```

## Testing and Validation

Mechanics implementations should be validated against:

### Analytical Tests
- **Two-cell collision** - Verify repulsion forces
- **Adhesion equilibrium** - Check steady-state separation
- **Energy conservation** - Spring oscillations

### Benchmark Problems
- **Random packing** - Dense cell configurations
- **Tissue spreading** - Adhesion-driven expansion
- **Compression** - Response to external forces

### Example: Two-Cell Collision

```cpp
TEST(Micromechanics, TwoCellRepulsion) {
    // Setup: Two overlapping cells
    auto solver = micromechanics_solver();
    solver.add_cell({0, 0, 0}, radius=10);
    solver.add_cell({15, 0, 0}, radius=10);  // Overlap of 5 units
    
    // Run mechanics
    solver.compute_forces();
    
    // Verify: Force magnitude
    real_t overlap = 5.0;
    real_t expected_force = k_rep * overlap;
    EXPECT_NEAR(solver.get_force(0).x, -expected_force, 1e-6);
    EXPECT_NEAR(solver.get_force(1).x, +expected_force, 1e-6);
}
```

## Future Directions

Planned enhancements to the mechanics module:

### Advanced Force Models
- **Viscoelastic cells** - Time-dependent deformation
- **Active forces** - Cell motility and migration
- **Substrate mechanics** - Extracellular matrix effects

### Numerical Methods
- **Implicit integration** - Larger timesteps
- **Contact mechanics** - Friction and sliding
- **Fluid coupling** - Stokes flow around cells

### Performance
- **GPU kernels** - Full GPU implementation
- **Adaptive timestepping** - Variable $\Delta t$
- **Domain decomposition** - Distributed parallelism

## Example Usage

```cpp
#include <mechanics/micromechanics/solver.h>

using namespace physicore::mechanics::micromechanics;

// Create mechanics solver
auto mechanics = micromechanics_solver();
mechanics.set_repulsion_coefficient(10.0);
mechanics.set_adhesion_coefficient(5.0);
mechanics.set_timestep(0.1);

// Add cells
for (int i = 0; i < 100; ++i) {
    real_t x = random(-100, 100);
    real_t y = random(-100, 100);
    real_t z = random(-100, 100);
    mechanics.add_cell({x, y, z}, radius=10.0);
}

// Simulation loop
for (int step = 0; step < 1000; ++step) {
    mechanics.run_single_timestep();
    
    if (step % 10 == 0) {
        mechanics.serialize_state(step * dt);
    }
}
```

## Next Steps

- **[Common Module](Architecture-Common.md)** - Agent containers and SoA patterns
- **[Reactions-Diffusion Module](Architecture-Diffusion.md)** - Substrate coupling
- **[Phenotype Module](Architecture-Phenotype.md)** - Integration and orchestration

---

**See also:**
- [Architecture Overview](Architecture.md)
- [Module Communication](Architecture-Phenotype.md#module-communication)
