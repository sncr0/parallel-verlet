// integrator/VerletIntegrator.h

#ifndef VERLET_INTEGRATOR_H
#define VERLET_INTEGRATOR_H

#include "../system/System.h"
#include <array>
#include <vector>

class VerletIntegrator {
public:
    explicit VerletIntegrator(double timestep);
    void step(System& system);

private:
    double timestep;
    
    // Parameters for Lennard-Jones potential
    const double epsilon = 1.0; // Depth of the potential well
    const double sigma = 1.0;   // Distance at which the potential is minimum
};

#endif
