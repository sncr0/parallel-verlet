// context/Context.cpp

#include "Context.h"
#include "../logging/Verbose.h"

Context::Context(System& sys, VerletIntegrator& integ) : system(sys), integrator(integ) {}

void Context::runSimulation(int steps) {
    for (int i = 0; i < steps; ++i) {
        integrator.step(system);
        if (steps % 10000 == 0) {
                for (size_t j = 0; j < system.getNumParticles(); ++j) {
                    auto pos = system.getParticle(j).getPosition();
                    auto vel = system.getParticle(j).getVelocity();
                    VERBOSE("\nTime step %d\n", i);
                    VERBOSE("Particle %zu position: (%f, %f, %f)\n", j, pos[0], pos[1], pos[2]);
                    VERBOSE("Particle %zu velocity: (%f, %f, %f)\n", j, vel[0], vel[1], vel[2]);
                }
        }
    }
}
