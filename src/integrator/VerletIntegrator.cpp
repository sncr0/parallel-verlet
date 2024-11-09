// integrator/VerletIntegrator.cpp

#include "VerletIntegrator.h"
#include <cmath> // for pow and sqrt

VerletIntegrator::VerletIntegrator(double timestep) : timestep(timestep) {}

void VerletIntegrator::step(System& system) {
    size_t numParticles = system.getNumParticles();
    
    // Temporary arrays for forces
    std::vector<std::array<double, 3>> forces(numParticles, {0.0, 0.0, 0.0});

    // Calculate Lennard-Jones forces between each pair of particles
    for (size_t i = 0; i < numParticles; ++i) {
        Particle& particle1 = system.getParticle(i);
        std::array<double, 3> pos1 = particle1.getPosition();

        for (size_t j = i + 1; j < numParticles; ++j) {
            Particle& particle2 = system.getParticle(j);
            std::array<double, 3> pos2 = particle2.getPosition();

            // Calculate the distance vector between particles i and j
            std::array<double, 3> r;
            double rSquared = 0.0;
            for (int k = 0; k < 3; ++k) {
                r[k] = pos2[k] - pos1[k];
                rSquared += r[k] * r[k];
            }
            double rMagnitude = std::sqrt(rSquared);

            // Avoid division by zero in case particles overlap
            if (rMagnitude == 0.0) continue;

            // Calculate Lennard-Jones force magnitude
            double sr6 = std::pow(sigma / rMagnitude, 6);
            double forceMagnitude = 24 * epsilon * (2 * sr6 * sr6 - sr6) / rSquared;

            // Update forces for particle i and j
            for (int k = 0; k < 3; ++k) {
                double force = forceMagnitude * r[k];
                forces[i][k] += force;     // Force on particle i due to j
                forces[j][k] -= force;     // Newton's Third Law
            }
        }
    }

    // Update positions using Verlet integration
    for (size_t i = 0; i < numParticles; ++i) {
        Particle& particle = system.getParticle(i);
        std::array<double, 3> pos = particle.getPosition();
        std::array<double, 3> vel = particle.getVelocity();

        // Update position: x(t+dt) = x(t) + v(t)*dt + (0.5 * F/m) * dt^2
        for (int k = 0; k < 3; ++k) {
            pos[k] += vel[k] * timestep + 0.5 * forces[i][k] * timestep * timestep / particle.getMass();
        }
        particle.setPosition(pos[0], pos[1], pos[2]);

        // Update velocity: v(t+dt) = v(t) + (0.5 * (F(t) + F(t+dt)) / m) * dt
        for (int k = 0; k < 3; ++k) {
            vel[k] += 0.5 * forces[i][k] * timestep / particle.getMass();
        }
        particle.setVelocity(vel[0], vel[1], vel[2]);
    }
}
