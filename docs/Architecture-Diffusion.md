---
layout: default
title: Reactions-Diffusion Module
parent: Architecture
nav_order: 2
has_children: true
description: "Substrate transport and reaction kinetics using reaction-diffusion PDEs in PhysiCore"
---

# Reactions-Diffusion Module
{: .no_toc }

The reactions-diffusion module handles substrate transport and reaction kinetics in multicellular simulations. It solves partial differential equations (PDEs) governing how chemical species diffuse through tissue, react with each other, and are produced or consumed by cells.

## Table of Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Overview

Substrate dynamics operate on timescales of **seconds to minutes**, making them one of the fastest processes in multicellular simulations. The reactions-diffusion module provides:

- **Diffusion** - Substrate transport via Brownian motion
- **Reactions** - Chemical transformations between species
- **Decay** - Natural degradation of substrates
- **Source/sink terms** - Cell secretion and uptake

**Namespace:** `physicore::reactions-diffusion`

**Location:** `reactions-diffusion/`

## Mathematical Foundation

### Governing Equation

The module solves the reaction-diffusion-decay equation for multiple substrates in 3D:

$$
\frac{\partial \rho_i}{\partial t} = D_i \nabla^2 \rho_i - \lambda_i \rho_i + S_i(\mathbf{x}, t)
$$

Where:
- $\rho_i(\mathbf{x}, t)$ - Concentration of substrate $i$ at position $\mathbf{x}$ and time $t$
- $D_i$ - Diffusion coefficient (substrate mobility)
- $\lambda_i$ - Decay rate (natural degradation)
- $S_i(\mathbf{x}, t)$ - Source/sink term (cell secretion/uptake)
- $\nabla^2$ - Laplacian operator (spatial second derivative)

### Physical Interpretation

**Diffusion term** ($D_i \nabla^2 \rho_i$):
- Models Brownian motion of molecules
- Higher diffusion coefficients → faster spread
- Driven by concentration gradients

**Decay term** ($-\lambda_i \rho_i$):
- Natural degradation or consumption
- Exponential decay: $\rho(t) = \rho_0 e^{-\lambda t}$
- Half-life: $t_{1/2} = \ln(2)/\lambda$

**Source/sink term** ($S_i(\mathbf{x}, t)$):
- Cell secretion (positive source)
- Cell uptake (negative sink)
- Spatially localized at cell positions

### Multi-Substrate Systems

For $N$ substrates, the system becomes:

$$
\frac{\partial \rho_i}{\partial t} = D_i \nabla^2 \rho_i - \lambda_i \rho_i + S_i(\mathbf{x}, t) + R_i(\rho_1, \ldots, \rho_N)
$$

Where $R_i(\rho_1, \ldots, \rho_N)$ represents chemical reactions between species.

## Module Interface Contract

All implementations of the reactions-diffusion module must implement the `timestep_executor` interface from the Common module:

```cpp
#include <common/timestep_executor.h>

namespace physicore::reactions_diffusion {

class diffusion_solver : public physicore::common::timestep_executor {
public:
    // Run one timestep of diffusion
    void run_single_timestep() override;
    
    // Serialize substrate concentrations
    void serialize_state(real_t current_time) override;
    
    // Module-specific interface
    virtual void set_diffusion_coefficient(std::size_t substrate_index, real_t D) = 0;
    virtual void set_decay_rate(std::size_t substrate_index, real_t lambda) = 0;
    virtual real_t get_concentration(std::size_t substrate_index, 
                                     real_t x, real_t y, real_t z) const = 0;
};

} // namespace physicore::reactions_diffusion
```

## Discretization Methods

The reactions-diffusion module can be implemented using various numerical methods:

### Finite Volume Method (FVM)
- Domain discretized into control volumes
- Concentrations stored at cell centers
- Fluxes computed at cell faces
- Conservative: mass is preserved
- **Implementation:** [BioFVM](BioFVM.md)

### Finite Difference Method (FDM)
- Regular grid discretization
- Direct approximation of derivatives
- Simple implementation
- **Status:** Future implementation

### Finite Element Method (FEM)
- Unstructured mesh support
- Complex geometries
- Higher-order accuracy
- **Status:** Future implementation

## Boundary Conditions

Implementations must support boundary conditions:

### Dirichlet Boundaries
Fix concentration at boundaries:
$$
\rho_i(\mathbf{x}_{\text{boundary}}, t) = \rho_i^{\text{boundary}}
$$

**Use cases:**
- Constant external concentration
- Infinite reservoirs
- Controlled environments

### Neumann Boundaries
Fix flux at boundaries:
$$
\nabla \rho_i \cdot \mathbf{n} = 0 \quad \text{(no-flux)}
$$

**Use cases:**
- Closed systems
- Isolated domains
- Zero-flux boundaries

## Time Integration

The module typically uses operator splitting for time integration:

1. **Diffusion step** - Solve $\frac{\partial \rho}{\partial t} = D \nabla^2 \rho$
2. **Reaction/decay step** - Solve $\frac{\partial \rho}{\partial t} = -\lambda \rho + S$

This allows different numerical methods for each operator, optimizing stability and accuracy.

## Coupling with Cells

The reactions-diffusion module interacts with cell agents through:

### Cell Secretion
Cells add substrates to their local voxels:
$$
S_i^{\text{secretion}}(\mathbf{x}, t) = \sum_{k} s_{i,k} \delta(\mathbf{x} - \mathbf{x}_k)
$$

Where:
- $s_{i,k}$ - Secretion rate of substrate $i$ by cell $k$
- $\mathbf{x}_k$ - Position of cell $k$
- $\delta(\mathbf{x} - \mathbf{x}_k)$ - Dirac delta (localized source)

### Cell Uptake
Cells consume substrates proportional to local concentration:
$$
S_i^{\text{uptake}}(\mathbf{x}, t) = -\sum_{k} u_{i,k} \rho_i(\mathbf{x}_k, t) \delta(\mathbf{x} - \mathbf{x}_k)
$$

Where:
- $u_{i,k}$ - Uptake rate of substrate $i$ by cell $k$

## Performance Considerations

Diffusion is often the computational bottleneck in multicellular simulations:

### Timestep Constraints
The diffusion equation imposes stability constraints:
$$
\Delta t \leq \frac{(\Delta x)^2}{6D}
$$

For typical values ($D \sim 10^5$ µm²/min, $\Delta x = 20$ µm):
$$
\Delta t \leq 0.67 \text{ min}
$$

### Computational Cost
- **Grid size:** $N_x \times N_y \times N_z$ voxels
- **Substrates:** $N_s$ chemical species
- **Operations per timestep:** $\mathcal{O}(N_x N_y N_z N_s)$

### Optimization Strategies
1. **Adaptive timestepping** - Larger steps when possible
2. **Sparse storage** - Only active regions
3. **Multigrid methods** - Faster convergence
4. **GPU acceleration** - Parallel voxel updates
5. **Operator splitting** - Optimize each sub-step

## Implementations

PhysiCore currently provides one production-ready implementation:

### BioFVM Implementation

**BioFVM** (Biological Finite Volume Method) is a highly optimized finite volume solver for reaction-diffusion PDEs.

**Key Features:**
- 3D Cartesian mesh discretization
- Multiple substrates with independent properties
- Dirichlet boundary conditions
- Efficient source/sink handling
- Pluggable solver backends (OpenMP, Thrust)

👉 **[Learn more about BioFVM](BioFVM.md)**

### Future Implementations

Planned implementations include:
- **Adaptive mesh refinement** - Higher resolution near cells
- **Unstructured meshes** - Complex tissue geometries
- **Spectral methods** - Higher-order accuracy

## Integration with Other Modules

The reactions-diffusion module integrates with:

### Mechanics Module
- **Input:** Cell positions and volumes from mechanics
- **Output:** Substrate concentrations influence cell behaviors
- **Coupling:** Source/sink terms localized at cell positions

### Phenotype Module
- **Input:** Cell phenotype determines secretion/uptake rates
- **Output:** Substrate levels trigger phenotype transitions
- **Coupling:** Bidirectional feedback between substrates and behaviors

See [Phenotype Module](Architecture-Phenotype.md#module-communication) for orchestration details.

## Example Usage

```cpp
#include <biofvm/microenvironment.h>
#include <biofvm/solver_registry.h>

using namespace physicore::reactions_diffusion::biofvm;

// Create microenvironment
auto microenv = microenvironment_builder()
    .set_domain(-500, 500, -500, 500, -500, 500)
    .set_voxel_size(20.0)
    .add_substrate("oxygen", 1e5, 0.1)     // D=1e5 µm²/min, λ=0.1/min
    .add_substrate("glucose", 8e4, 0.05)
    .build();

// Select solver backend
auto solver = solver_registry::instance().get("openmp");
microenv.set_solver(std::move(solver));

// Simulation loop
for (int step = 0; step < 1000; ++step) {
    microenv.run_single_timestep();
    
    if (step % 10 == 0) {
        microenv.serialize_state(step * dt);
    }
}
```

## Testing

Implementations should be validated against:
- **Analytical solutions** - Simple geometries and boundary conditions
- **Benchmark problems** - Standard test cases from literature
- **Conservation laws** - Mass conservation in closed systems
- **Convergence rates** - Spatial and temporal accuracy

## Next Steps

- **[BioFVM Implementation](BioFVM.md)** - Detailed guide to the finite volume solver
- **[Common Module](Architecture-Common.md)** - Understanding the `timestep_executor` interface
- **[Phenotype Module](Architecture-Phenotype.md)** - How diffusion integrates with complete simulations

---

**See also:**
- [Architecture Overview](Architecture.md)
- [Kernel Implementations](Architecture.md#kernel-implementations)
