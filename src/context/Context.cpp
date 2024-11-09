// context/Context.cpp

#include "Context.h"

Context::Context(System& sys, VerletIntegrator& integ) : system(sys), integrator(integ) {}

void Context::runSimulation(int steps) {
    for (int i = 0; i < steps; ++i) {
        integrator.step(system);
    }
}
