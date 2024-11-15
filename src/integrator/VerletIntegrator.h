// integrator/VerletIntegrator.h

#ifndef VERLET_INTEGRATOR_H
#define VERLET_INTEGRATOR_H

#include "../system/System.h"
#include <array>
#include <vector>
#include <memory>
#include "../forces/Force.h"

class VerletIntegrator {
public:
    explicit VerletIntegrator(double timestep);
    void step(System& system);
    void addForce(std::shared_ptr<Force> force);


private:
    double timestep;
    
    // Parameters for Lennard-Jones potential
    const double epsilon = 1.0; // Depth of the potential well
    const double sigma = 1.0;   // Distance at which the potential is minimum

    std::vector<std::shared_ptr<Force>> forces; // Collection of forces

};

#endif
