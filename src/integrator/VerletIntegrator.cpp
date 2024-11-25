// integrator/VerletIntegrator.cpp

#include "VerletIntegrator.h"
#include <cmath> // for pow and sqrt
#include "../logging/Verbose.h"
#include "../logging/Chronometer.h"
#include <omp.h>


/* Velocity Verlet Integrator for Molecular Dynamics
 *
 * Physical Description:
 * This integrator evolves the system forward in time by updating particle positions and velocities
 * based on the forces acting on them. It uses the Velocity Verlet algorithm, which is symplectic
 * (meaning it preserves phase-space volume) and time-reversible.
 *
 * The algorithm splits each timestep into multiple phases to achieve second-order accuracy:
 * 1. Position update using current velocity and acceleration (from forces):
 *    r(t + Δt) = r(t) + v(t)Δt + (1/2)a(t)(Δt)²
 *    This comes from Taylor expansion of position around t
 *
 * 2. Velocity update in two half-steps:
 *    First half:  v(t + Δt/2) = v(t) + (1/2)a(t)Δt
 *    Second half: v(t + Δt) = v(t + Δt/2) + (1/2)a(t + Δt)Δt
 *    The split update allows using both old and new forces for better accuracy
 *
 * Key Physical Properties:
 * - Symplectic: Preserves phase-space volume, leading to good energy conservation
 * - Time-reversible: Running simulation backwards gives same trajectory
 * - Second-order accurate: Local error O(Δt³) in positions, O(Δt²) in velocities
 * - Synchronized: Positions and velocities are computed at the same time points
 *
 * Force Handling:
 * - Forces are stored between timesteps to avoid double calculation
 * - F = ma is used to get acceleration: a = F/m
 * - Forces must be conservative (derived from potential energy) for energy conservation
 *
 * Implementation Details:
 * - Uses member variable old_forces to store forces between steps
 * - first_step flag ensures proper initialization
 * - Handles multiple force components through the forces vector
 */

VerletIntegrator::VerletIntegrator(double ts, 
                                    ThreadManager& tm,
                                    Chronometer& chrono) : timestep(ts), 
                                                            thread_manager(tm),
                                                            chronometer(chrono) {}

void VerletIntegrator::addForce(std::shared_ptr<Force> force) {
    forces.push_back(force);
}

void VerletIntegrator::step(System& system) {
    chronometer.start("step");
    size_t numParticles = system.getNumParticles();

    // Initialize old_forces if this is the first step
    if (first_step) {
        old_forces.resize(numParticles, {0.0, 0.0, 0.0});
        // Calculate initial forces
        for (const auto& force : forces) {
            if (force->num_threads == 0) {
                force->compute(system, old_forces, thread_manager, chronometer);
            } else {
                force->compute_parallel(system, old_forces, thread_manager, chronometer);
            }
        }
        first_step = false;
    }
   
    // Step 1: Update positions and half-step velocities using old forces
    for (size_t i = 0; i < numParticles; ++i) {
        Particle& particle = system.getParticle(i);
        auto pos = particle.getPosition();
        auto vel = particle.getVelocity();
        auto mass = particle.getMass();

        for (int k = 0; k < 3; ++k) {
            // Update position using current velocity and old forces
            pos[k] += vel[k] * timestep + 0.5 * old_forces[i][k] * timestep * timestep / mass;
            // First half of velocity update using old forces
            vel[k] += 0.5 * old_forces[i][k] * timestep / mass;
        }
        particle.setPosition(pos[0], pos[1], pos[2]);
        particle.setVelocity(vel[0], vel[1], vel[2]);
    }

    // Step 2: Calculate new forces at new positions
    std::vector<std::array<double, 3>> new_forces(numParticles, {0.0, 0.0, 0.0});
    #pragma omp parallel num_threads(2)
    {
        std::vector<std::array<double, 3>> new_forces_thread(numParticles, {0.0, 0.0, 0.0});
        #pragma omp for
        for (const auto& force : forces) {
            if (force->num_threads == 0) {
                force->compute(system, new_forces_thread, thread_manager, chronometer);
            } else {
                force->compute_parallel(system, new_forces_thread, thread_manager, chronometer);
            }
            // force->compute(system, new_forces_thread, thread_manager, chronometer);
        }

        #pragma omp critical
        {
            for (size_t i = 0; i < numParticles; ++i) {
                for (int k = 0; k < 3; ++k) {
                    new_forces[i][k] += new_forces_thread[i][k];
                }
            }
        }
    }

    // Step 3: Complete velocity update with new forces
    for (size_t i = 0; i < numParticles; ++i) {
        Particle& particle = system.getParticle(i);
        auto vel = particle.getVelocity();
        auto mass = particle.getMass();

        for (int k = 0; k < 3; ++k) {
            // Second half of velocity update using new forces
            vel[k] += 0.5 * new_forces[i][k] * timestep / mass;
        }
        particle.setVelocity(vel[0], vel[1], vel[2]);
    }

    // Store new forces for next timestep
    old_forces = new_forces;

    chronometer.end("step");
    chronometer.printTiming("step", "ms", 0);
}
