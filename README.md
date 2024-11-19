# Project Overview: Molecular Dynamics Simulation Framework

This project implements a molecular dynamics simulation framework, including components for managing the simulation context, applying physical forces, integrating equations of motion, and outputting results.

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
