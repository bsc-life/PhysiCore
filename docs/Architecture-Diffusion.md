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

## Implementations

PhysiCore currently provides one production-ready implementation:

### BioFVM Implementation

**BioFVM** (Biological Finite Volume Method) is a highly optimized finite volume solver for reaction-diffusion PDEs.

**Key Features:**
- 3D Cartesian mesh discretization
- Multiple substrates with independent properties
- Dirichlet boundary conditions
- Efficient source/sink handling

👉 **[Learn more about BioFVM](BioFVM.md)**

<!-- ## Integration with Other Modules

The reactions-diffusion module integrates with:

### Mechanics Module
- **Input:** Cell positions and volumes from mechanics
- **Output:** Substrate concentrations influence cell behaviors
- **Coupling:** Source/sink terms localized at cell positions

### Phenotype Module
- **Input:** Cell phenotype determines secretion/uptake rates
- **Output:** Substrate levels trigger phenotype transitions
- **Coupling:** Bidirectional feedback between substrates and behaviors

See [Phenotype Module](Architecture-Phenotype.md#module-communication) for orchestration details. -->

## Next Steps

- **[BioFVM Implementation](BioFVM.md)** - Detailed guide to the finite volume solver
- **[Common Module](Architecture-Common.md)** - Understanding the `timestep_executor` interface
<!-- - **[Phenotype Module](Architecture-Phenotype.md)** - How diffusion integrates with complete simulations -->

---

**See also:**
- [Architecture Overview](Architecture.md)
