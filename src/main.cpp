// main.cpp

#include "system/System.h"
// #include "context/Context.h"
// #include "integrator/VerletIntegrator.h"
#include <iostream>

int main() {
    System system;
    system.addParticle(1.0, 0.0, 0.0, 0.0);
    system.addParticle(1.0, 1.0, 1.0, 1.0);

    // VerletIntegrator integrator(0.01);
    // Context context(system, integrator);

    // // Run the simulation
    // context.runSimulation(10);

    // Output the final positions of particles
    for (size_t i = 0; i < system.getNumParticles(); ++i) {
        auto pos = system.getParticle(i).getPosition();
        std::cout << "Particle " << i << " final position: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
    }

    return 0;
}
