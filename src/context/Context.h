// context/Context.h

#ifndef CONTEXT_H
#define CONTEXT_H

#include "../system/System.h"
#include "../integrator/VerletIntegrator.h"

class Context {
public:
    Context(System& system, VerletIntegrator& integrator);
    void runSimulation(int steps);

private:
    System& system;
    VerletIntegrator& integrator;
};

#endif
