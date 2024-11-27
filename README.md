# Project Overview: Molecular Dynamics Simulation Framework

This project implements a molecular dynamics simulation framework, including components for managing the simulation context, applying physical forces, integrating equations of motion, and outputting results.

### **Compiling the Project with Make**

The project includes a `Makefile` to simplify the compilation process. You can use the `make` command to build the executables and manage dependencies.

```bash
make
```


## Performance Benchmarking

### Available Benchmark Scripts

The project provides several benchmark scripts for performance analysis:

- `benchmark_bonds.sh`
- `benchmark_dispersion.sh`
- `benchmark_electrostatic.sh`

### Running Benchmarks

```bash
# General benchmark execution format
./benchmark_<type> -name
```

### Benchmark Output

- Results are saved to `results/data`
- Performance visualizations available in `results/plot_speedups.ipynb`

## Executable Options

### Electrostatic Interactions Executable

#### Parallel Execution
```bash
./bin/benchmark_electrostatic -n <system_size> -s <steps> -v <verbose> -w <write> -e <threads>
```

#### Sequential Execution
```bash
./bin/benchmark_electrostatic -n <system_size> -s <steps> -v <verbose> -w <write>
./bin/benchmark_electrostatic -n <system_size> -s <steps> -v <verbose> -w <write> -e 0
```

### Dispersion Interactions Executable

#### Parallel Execution
```bash
./bin/benchmark_dispersion -n <system_size> -s <steps> -v <verbose> -w <write> -d <threads>
```

#### Sequential Execution
```bash
./bin/benchmark_dispersion -n <system_size> -s <steps> -v <verbose> -w <write>
./bin/benchmark_dispersion -n <system_size> -s <steps> -v <verbose> -w <write> -d 0
```

### Bonded Interactions Executable

#### Parallel Execution
```bash
./bin/benchmark_bonds -n <system_size> -s <steps> -v <verbose> -w <write> -h <threads>
```

#### Sequential Execution
```bash
./bin/benchmark_bonds n <system_size> -s <steps> -v <verbose> -w <write>
./bin/benchmark_bonds -n <system_size> -s <steps> -v <verbose> -w <write> -h 0

```

## Command-Line Option Details

| Option | Description |
|--------|-------------|
| `-n`   | System size (number of particles) |
| `-s`   | Number of simulation steps |
| `-v`   | Verbose output mode |
| `-w`   | Write output to file |
| `-e`   | Threads for electrostatic calculations |
| `-d`   | Threads for dispersion calculations |
| `-h`   | Threads for bonded interactions |

## Example Usage

### Parallel Electrostatic Simulation
```bash
# Run with 1000 particles, 500 steps, 4 threads
./bin/benchmark_electrostatic -n 1000 -s 500 -e 4
```

### Sequential Electrostatic Simulation
```bash
# Run with 1000 particles, 500 steps
./bin/benchmark_electrostatic -n 1000 -s 500
```


## **Directory Structure and Functional Overview**

### **Main File**
#### `main.cpp`
- **Purpose**: Entry point for the program.
- **Responsibilities**:
  - Initializes the **`System`** by defining particles and their properties.
  - Configures the **`Integrator`** (e.g., `VerletIntegrator`) and registers physical **`Forces`**.
  - Sets up the simulation **`Context`**, linking the `System`, `Integrator`, and forces.
  - Runs the simulation for a specified number of steps.
  - Outputs the simulation results, such as particle positions, to the console or a file.

---

### **Context**
#### `context/Context.h` and `context/Context.cpp`
- **Purpose**: Serves as the simulation controller, managing interactions between the `System`, `Forces`, and `Integrator`.
- **Responsibilities**:
  - Maintains the current simulation state.
  - Manages the integration process step-by-step.
  - 
---

### **Forces**
#### `forces/Force.h` and `forces/Force.cpp`
- **Purpose**: Abstract base class for defining physical forces acting on particles.
- **Responsibilities**:
  - Provides a generic interface for all forces (e.g., `applyForce`).
  - Allows for custom force implementations by extending this base class.

#### `forces/LennardJonesForce.h` and `forces/LennardJonesForce.cpp`
- **Purpose**: Implements the Lennard-Jones potential for non-bonded particle interactions.
- **Responsibilities**:
  - Calculates attractive and repulsive forces between particle pairs.

#### `forces/HarmonicBondForce.h` and `forces/HarmonicBondForce.cpp`
- **Purpose**: Implements a harmonic potential to simulate bonded particle interactions.
- **Responsibilities**:
  - Computes forces and potential energy for bonded particle pairs.

---

### **Integrator**
#### `integrator/VerletIntegrator.h` and `integrator/VerletIntegrator.cpp`
- **Purpose**: Implements the Verlet integration algorithm for time-stepping.
- **Responsibilities**:
  - Updates particle positions and velocities based on forces.
  - Ensures numerical stability during the simulation.

---

### **System**
#### `system/System.h` and `system/System.cpp`
- **Purpose**: Represents the collection of particles in the simulation.
- **Responsibilities**:
  - Manages particle properties (e.g., mass, position, velocity).
  - Provides methods to add, access, and update particles.

#### `system/Particle.h` and `system/Particle.cpp`
- **Purpose**: Represents an individual particle in the system.
- **Responsibilities**:
  - Encapsulates particle attributes, such as position, velocity, and mass.
  - Provides methods to retrieve or update particle properties.

---

### **IO**
#### `io/XYZWriter.h` and `io/XYZWriter.cpp`
- **Purpose**: Handles writing simulation trajectories to files in the XYZ format.
- **Responsibilities**:
  - Outputs particle positions for each simulation frame.
  - Ensures compatibility with visualization tools like VMD.

---

### **Logging**
#### `logging/Verbose.h` and `logging/Verbose.cpp`
- **Purpose**: Provides utilities for managing verbose logging output.
- **Responsibilities**:
  - Allows enabling or disabling detailed log messages.
  - Outputs simulation information for debugging and monitoring purposes.


