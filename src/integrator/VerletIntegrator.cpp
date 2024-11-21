// integrator/VerletIntegrator.cpp

#include "VerletIntegrator.h"
#include <cmath> // for pow and sqrt
#include "../logging/Verbose.h"

VerletIntegrator::VerletIntegrator(double ts, 
                                    ThreadManager& tm) : timestep(ts), 
                                                            thread_manager(tm) {}

void VerletIntegrator::addForce(std::shared_ptr<Force> force) {
    forces.push_back(force);
}

void VerletIntegrator::step(System& system) {
    size_t numParticles = system.getNumParticles();

    // for (size_t i = 0; i < system.getNumParticles(); ++i) {
    //     auto pos = system.getParticle(i).getPosition();
    //     std::cout << "Particle " << i << " position: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
    // }
    
    // Temporary arrays for forces
    std::vector<std::array<double, 3>> forcesArray(numParticles, {0.0, 0.0, 0.0});

    for (const auto& force : forces) {
        force->compute(system, forcesArray, thread_manager);
    }

    // Integrate positions and velocities
    for (size_t i = 0; i < numParticles; ++i) {
        Particle& particle = system.getParticle(i);
        auto pos = particle.getPosition();
        auto vel = particle.getVelocity();
        auto mass = particle.getMass();

        for (int k = 0; k < 3; ++k) {
            pos[k] += vel[k] * timestep + 0.5 * forcesArray[i][k] * timestep * timestep / mass;
            vel[k] += 0.5 * forcesArray[i][k] * timestep / mass;
        }
        particle.setPosition(pos[0], pos[1], pos[2]);
        particle.setVelocity(vel[0], vel[1], vel[2]);
    }
}
