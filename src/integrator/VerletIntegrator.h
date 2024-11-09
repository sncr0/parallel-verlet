// integrator/VerletIntegrator.h

#ifndef VERLET_INTEGRATOR_H
#define VERLET_INTEGRATOR_H

#include "../system/System.h"

class VerletIntegrator {
public:
    explicit VerletIntegrator(double timestep);
    void step(System& system);

private:
    double timestep;
};

#endif
