---
layout: default
title: Phenotype Module
parent: Architecture
nav_order: 4
description: "Integration layer coordinating diffusion and mechanics modules with biological phenotype behaviors"
---

# Phenotype Module
{: .no_toc }

The phenotype module is PhysiCore's integration layer that coordinates the reactions-diffusion and mechanics modules into complete multicellular simulations. It defines cell phenotype behaviors, orchestrates timesteps across modules, and manages simulation state.

## Table of Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Overview

The phenotype module operates at the **biological behavior timescale** (minutes to hours), coordinating faster substrate diffusion with slower mechanical interactions. It provides:

- **Module orchestration** - Coordinates diffusion and mechanics timesteps
- **Phenotype behaviors** - Defines cell state transitions and rules
- **Biological coupling** - Links substrate levels to cell behaviors
- **State management** - Serialization and checkpointing
- **Simulation wiring** - Main loop and integration logic

**Namespace:** `physicore::phenotype`

**Location:** `phenotype/`

## Architectural Role

The phenotype module sits at the top of the module hierarchy:

```
┌─────────────────────────────────────────┐
│           Phenotype Layer               │
│  (Coordinates simulation, biology)      │
└───────────────┬─────────────────────────┘
                │
       ┌────────┴────────┐
       │                 │
┌──────▼──────┐   ┌──────▼──────┐
│  Diffusion  │   │  Mechanics  │
│   Module    │◄──┤   Module    │
└─────────────┘   └─────────────┘
  Concentrations    Cell Positions
```

**Responsibilities:**
1. **Coordinate** diffusion and mechanics modules
2. **Define** cell phenotype states and transitions
3. **Couple** substrate concentrations to cell behaviors
4. **Manage** simulation time and output
5. **Wire** all components into executable simulations

## Module Communication

The phenotype module orchestrates bidirectional communication between modules:

### Diffusion → Phenotype → Mechanics

**Flow:**
1. Diffusion module computes substrate concentrations
2. Phenotype reads concentrations at cell positions
3. Substrate levels influence cell behaviors:
   - Secretion/uptake rates
   - Growth rates
   - Migration speed
   - Phenotype transitions
4. Updated behaviors affect mechanics forces

**Example:**
```cpp
// Read oxygen concentration at cell position
real_t oxygen = diffusion.get_concentration("oxygen", cell.position());

// Oxygen affects cell growth
if (oxygen < threshold) {
    cell.set_growth_rate(0.0);  // Hypoxic, no growth
} else {
    cell.set_growth_rate(0.05);  // Normal growth
}
```

### Mechanics → Phenotype → Diffusion

**Flow:**
1. Mechanics module updates cell positions
2. Phenotype determines secretion/uptake based on position
3. Diffusion module uses cell positions for source/sink terms

**Example:**
```cpp
// After mechanics update
mechanics.run_single_timestep();

// Update diffusion source terms based on new positions
for (auto& cell : cells) {
    diffusion.add_source("oxygen", 
        cell.position(), 
        cell.secretion_rate());
}
```

### Phenotype Feedback Loops

Substrate levels trigger cell state changes:

```cpp
class CellPhenotype {
    enum State { Normoxic, Hypoxic, Necrotic };
    State state = Normoxic;
    
    void update(real_t oxygen) {
        switch (state) {
            case Normoxic:
                if (oxygen < hypoxia_threshold) {
                    state = Hypoxic;
                    secretion_rate *= 0.5;  // Reduced metabolism
                }
                break;
            case Hypoxic:
                if (oxygen < necrotic_threshold) {
                    state = Necrotic;
                    secretion_rate = 0.0;  // Cell death
                }
                break;
        }
    }
};
```

## Timestep Orchestration

Different modules operate at different timescales. The phenotype module coordinates them:

### Timescale Hierarchy

| Module | Timescale | Typical $\Delta t$ |
|--------|-----------|-------------------|
| Diffusion | Seconds to minutes | 0.01 - 0.1 min |
| Mechanics | Minutes | 0.1 - 1.0 min |
| Phenotype | Hours | 1.0 - 10.0 min |

### Operator Splitting

The phenotype module uses **operator splitting** to integrate modules with different timesteps:

```cpp
class PhenotypeSolver : public physicore::common::timestep_executor {
    DiffusionSolver& diffusion;
    MechanicsSolver& mechanics;
    real_t phenotype_dt = 1.0;  // Phenotype timestep
    
public:
    void run_single_timestep() override {
        // 1. Run multiple diffusion substeps
        int diffusion_substeps = phenotype_dt / diffusion.get_dt();
        for (int i = 0; i < diffusion_substeps; ++i) {
            update_diffusion_sources();  // Cell secretion/uptake
            diffusion.run_single_timestep();
        }
        
        // 2. Run mechanics at phenotype timescale
        update_cell_forces();
        mechanics.run_single_timestep();
        
        // 3. Update cell phenotypes based on environment
        for (auto& cell : cells) {
            real_t oxygen = diffusion.get_concentration("oxygen", cell.position());
            cell.update_phenotype(oxygen);
        }
        
        current_time += phenotype_dt;
    }
};
```

### Adaptive Timestepping

Advanced implementations can adapt timesteps based on activity:

```cpp
real_t compute_adaptive_timestep() {
    // Limit based on cell movement
    real_t max_velocity = get_max_cell_velocity();
    real_t dt_mechanics = 0.1 * voxel_size / max_velocity;
    
    // Limit based on concentration gradients
    real_t max_gradient = diffusion.get_max_gradient();
    real_t dt_diffusion = 0.1 * voxel_size / (D * max_gradient);
    
    // Use minimum to ensure stability
    return std::min({dt_mechanics, dt_diffusion, dt_max});
}
```

## Phenotype Behavior Models

The phenotype module defines biological rules governing cell behavior:

### State Machine Approach

Cells transition between discrete states:

```cpp
class CellPhenotype {
    enum State {
        Quiescent,      // Resting, not dividing
        Proliferating,  // Actively growing
        Apoptotic,      // Programmed death
        Necrotic        // Dead
    };
    
    State state;
    real_t state_timer;
    
    void update(real_t oxygen, real_t nutrients) {
        switch (state) {
            case Quiescent:
                if (oxygen > proliferation_threshold && 
                    nutrients > proliferation_threshold) {
                    state = Proliferating;
                    state_timer = 0.0;
                }
                break;
                
            case Proliferating:
                state_timer += dt;
                if (state_timer > division_time) {
                    divide();
                    state = Quiescent;
                }
                if (oxygen < hypoxia_threshold) {
                    state = Quiescent;
                }
                break;
                
            case Apoptotic:
                state_timer += dt;
                if (state_timer > apoptosis_duration) {
                    remove_cell();
                }
                break;
        }
    }
};
```

### Continuous Rate Models

Alternative: Use continuous rates instead of discrete states:

```cpp
class CellPhenotype {
    real_t proliferation_rate;
    real_t apoptosis_rate;
    real_t necrosis_rate;
    
    void update(real_t oxygen) {
        // Rates depend on microenvironment
        proliferation_rate = hill_function(oxygen, K_prolif, n);
        necrosis_rate = (oxygen < necrotic_threshold) ? 1.0 : 0.0;
    }
    
    void stochastic_update(real_t dt) {
        // Stochastic events based on rates
        real_t p_divide = 1.0 - exp(-proliferation_rate * dt);
        real_t p_die = 1.0 - exp(-necrosis_rate * dt);
        
        if (random() < p_divide) divide();
        if (random() < p_die) die();
    }
};
```

## Module Interface Contract

The phenotype module implements the `timestep_executor` interface:

```cpp
#include <common/timestep_executor.h>

namespace physicore::phenotype {

class phenotype_solver : public physicore::common::timestep_executor {
public:
    // Run one phenotype timestep (coordinates all modules)
    void run_single_timestep() override;
    
    // Serialize complete simulation state
    void serialize_state(real_t current_time) override;
    
    // Module-specific interface
    virtual void register_diffusion_module(DiffusionSolver& diffusion) = 0;
    virtual void register_mechanics_module(MechanicsSolver& mechanics) = 0;
    virtual void add_cell_type(const CellTypeDefinition& type) = 0;
};

} // namespace physicore::phenotype
```

## Implementations

### PhysiCore Phenotype Implementation

**PhysiCore Phenotype** is the reference implementation of the phenotype module.

**Namespace:** `physicore::phenotype::physicore`

**Location:** `phenotype/physicore/`

**Key Features:**
- Coordinates BioFVM diffusion and micromechanics
- State machine-based phenotype behaviors
- Stochastic cell division and death
- VTK output for cells and substrates
- Configuration via XML

**Status:** Under active development

## Simulation Wiring

The phenotype module provides the main simulation executable:

### Basic Simulation Loop

```cpp
#include <phenotype/physicore/simulation.h>
#include <biofvm/microenvironment.h>
#include <mechanics/micromechanics/solver.h>

int main() {
    // 1. Create modules
    auto diffusion = create_biofvm_solver();
    auto mechanics = create_micromechanics_solver();
    
    // 2. Create phenotype coordinator
    auto phenotype = create_phenotype_solver();
    phenotype.register_diffusion_module(diffusion);
    phenotype.register_mechanics_module(mechanics);
    
    // 3. Define cell types
    CellTypeDefinition cancer_cell = {
        .name = "cancer",
        .division_time = 24.0 * 60.0,  // 24 hours
        .apoptosis_rate = 0.0,
        .oxygen_uptake = 10.0,
        .oxygen_secretion = 0.0
    };
    phenotype.add_cell_type(cancer_cell);
    
    // 4. Initialize cells
    for (int i = 0; i < 100; ++i) {
        phenotype.add_cell("cancer", random_position());
    }
    
    // 5. Simulation loop
    real_t max_time = 7 * 24 * 60;  // 7 days
    real_t dt = 1.0;  // 1 minute
    
    for (real_t t = 0; t < max_time; t += dt) {
        phenotype.run_single_timestep();
        
        if (fmod(t, 60.0) < dt) {  // Every hour
            phenotype.serialize_state(t);
        }
    }
    
    return 0;
}
```

### Configuration-Driven Simulations

Production simulations use XML configuration:

```xml
<!-- simulation_config.xml -->
<simulation>
    <domain>
        <x_min>-500</x_min>
        <x_max>500</x_max>
        <!-- ... -->
    </domain>
    
    <substrates>
        <substrate name="oxygen">
            <diffusion_coefficient>1e5</diffusion_coefficient>
            <decay_rate>0.1</decay_rate>
        </substrate>
    </substrates>
    
    <cell_types>
        <cell_type name="cancer">
            <division_time>1440</division_time>
            <oxygen_uptake>10.0</oxygen_uptake>
        </cell_type>
    </cell_types>
    
    <initial_conditions>
        <cell_population type="cancer" count="100">
            <distribution>random_ball</distribution>
            <radius>100</radius>
        </cell_population>
    </initial_conditions>
</simulation>
```

Load configuration:
```cpp
auto config = load_config("simulation_config.xml");
auto simulation = create_simulation(config);
simulation.run();
```

## State Serialization

The phenotype module coordinates output from all modules:

### VTK Output Format

```cpp
void phenotype_solver::serialize_state(real_t current_time) {
    // 1. Serialize cell data
    write_cell_vtk(current_time, cells);
    
    // 2. Serialize substrate data
    diffusion.serialize_state(current_time);
    
    // 3. Serialize metadata
    write_metadata(current_time, {
        {"num_cells", cells.size()},
        {"simulation_time", current_time}
    });
}
```

Output structure:
```
output/
├── cells_00000.vtu        # Cell positions, radii, types
├── cells_00001.vtu
├── substrate_00000.vtu    # Oxygen, glucose concentrations
├── substrate_00001.vtu
└── metadata.json
```

### Checkpoint/Restart

For long simulations, save complete state:

```cpp
void phenotype_solver::save_checkpoint(const std::string& filename) {
    Archive archive(filename);
    
    archive.write("time", current_time);
    archive.write("cells", cells);
    archive.write("diffusion", diffusion.get_state());
    archive.write("mechanics", mechanics.get_state());
}

void phenotype_solver::load_checkpoint(const std::string& filename) {
    Archive archive(filename);
    
    current_time = archive.read<real_t>("time");
    cells = archive.read<CellContainer>("cells");
    diffusion.set_state(archive.read("diffusion"));
    mechanics.set_state(archive.read("mechanics"));
}
```

## Performance Considerations

### Module Coupling Overhead

Minimize data transfer between modules:

```cpp
// Inefficient: Query diffusion for each cell individually
for (auto& cell : cells) {
    real_t oxygen = diffusion.get_concentration("oxygen", cell.position());
    cell.update(oxygen);
}

// Efficient: Batch query
auto oxygen_values = diffusion.get_concentrations_at_positions(cell_positions);
for (size_t i = 0; i < cells.size(); ++i) {
    cells[i].update(oxygen_values[i]);
}
```

### Timestep Ratios

Optimize the ratio of module timesteps:

```cpp
// Good: Power-of-2 ratios enable clean nesting
dt_phenotype = 1.0;
dt_mechanics = 0.5;    // 2× finer
dt_diffusion = 0.25;   // 4× finer

// Bad: Awkward ratios cause synchronization issues
dt_phenotype = 1.0;
dt_mechanics = 0.3;    // 3.33 substeps needed
dt_diffusion = 0.07;   // 14.29 substeps needed
```

## Testing and Validation

### Integration Tests

Test module interactions:

```cpp
TEST(PhenotypeIntegration, DiffusionMechanicsCoupling) {
    // Setup coupled simulation
    auto sim = create_test_simulation();
    sim.add_cell({0, 0, 0});
    
    // Run simulation
    sim.run_single_timestep();
    
    // Verify coupling
    auto cell_pos = sim.get_cell_position(0);
    auto oxygen = sim.get_oxygen_at(cell_pos);
    
    EXPECT_GT(oxygen, 0.0);  // Diffusion provided oxygen
    EXPECT_TRUE(sim.cell_moved(0));  // Mechanics updated position
}
```

### Biological Validation

Compare against experimental data:
- **Tumor growth curves** - Match spheroid expansion rates
- **Spatial patterns** - Verify necrotic core formation
- **Cell counts** - Proliferation/death balance

## Future Directions

Planned enhancements:

### Advanced Phenotype Models
- **Gene regulatory networks** - Intracellular signaling
- **Cell-cell communication** - Juxtacrine/paracrine signaling
- **Heterogeneous populations** - Multiple interacting cell types

### Optimization
- **GPU acceleration** - Parallelize phenotype updates
- **Load balancing** - Distribute cells across processes
- **Hierarchical timestepping** - Adaptive multi-rate integration

### Usability
- **Python bindings** - High-level scripting interface
- **Parameter calibration** - Automated fitting to data
- **Uncertainty quantification** - Sensitivity analysis

## Example: Complete Tumor Simulation

```cpp
#include <physicore/simulation.h>

int main() {
    // Configuration
    SimulationConfig config;
    config.domain = {-500, 500, -500, 500, -500, 500};
    config.max_time = 7 * 24 * 60;  // 7 days
    config.phenotype_dt = 1.0;
    
    // Create simulation
    auto sim = Simulation(config);
    
    // Add oxygen substrate
    sim.add_substrate("oxygen")
        .set_diffusion_coefficient(1e5)
        .set_decay_rate(0.1)
        .set_boundary_value(21.0);  // 21% atmospheric
    
    // Define cancer cell type
    sim.add_cell_type("cancer")
        .set_division_time(24 * 60)
        .set_apoptosis_rate(0.0)
        .set_oxygen_uptake(10.0)
        .set_hypoxia_threshold(5.0);
    
    // Initialize tumor
    sim.add_cell_sphere("cancer", {0, 0, 0}, radius=50, count=20);
    
    // Run simulation
    sim.run();
    
    // Analysis
    std::cout << "Final cell count: " << sim.get_num_cells() << "\n";
    std::cout << "Tumor radius: " << sim.compute_tumor_radius() << " μm\n";
    
    return 0;
}
```

## Next Steps

- **[Common Module](Architecture-Common.md)** - Understanding the `timestep_executor` interface
- **[Reactions-Diffusion Module](Architecture-Diffusion.md)** - Substrate transport details
- **[Mechanics Module](Architecture-Mechanics.md)** - Force calculations and position updates

---

**See also:**
- [Architecture Overview](Architecture.md)
- [BioFVM Implementation](BioFVM.md)
- [Installation Guide](Installation.md)
