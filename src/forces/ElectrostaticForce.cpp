#include "ElectrostaticForce.h"
#include <cmath> // for sqrt and pow
#include <omp.h>  // for OpenMP
#include "../logging/Verbose.h"
#include <iostream>
#include <chrono>
#include "../thread_manager/ThreadManager.h"

ElectrostaticForce::ElectrostaticForce(double coulombConstant)
    : coulombConstant(coulombConstant) {}

void ElectrostaticForce::addBond(size_t particle1Index, size_t particle2Index) {
    bonds.emplace_back(particle1Index, particle2Index);
}

void ElectrostaticForce::compute_parallel(System& system,
                                  std::vector<std::array<double, 3>>& forces,
                                  ThreadManager& thread_manager,
                                  Chronometer& chronometer) const {
    VERBOSE("Computing electrostatic forces\n");
    printf("Computing electrostatic forces in parallel\n");
    chronometer.start("electrostatic_force_glob");

    // Create thread-local forces arrays
    int num_threads = thread_manager.electrostatic_threads;
    size_t num_bonds = bonds.size();
    size_t num_particles = forces.size();

    std::vector<std::vector<std::array<double, 3>>> local_forces(num_threads,
        std::vector<std::array<double, 3>>(forces.size(), {0.0, 0.0, 0.0}));

    #pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();
        chronometer.start("electrostatic_force", thread_id);
        std::vector<double> thread_local_forces_x(num_particles, 0.0);
        std::vector<double> thread_local_forces_y(num_particles, 0.0);
        std::vector<double> thread_local_forces_z(num_particles, 0.0);
        // std::vector<std::array<double, 3>> thread_local_forces(num_particles, {0.0, 0.0, 0.0});
        
        #pragma omp for collapse(2)
        for (size_t particle1_index = 0; particle1_index < num_particles; ++particle1_index) {
            for (size_t particle2_index = particle1_index + 1; particle2_index < num_particles; ++particle2_index) {
                Particle& particle1 = system.getParticle(particle1_index);
                Particle& particle2 = system.getParticle(particle2_index);

                const std::array<double, 3>& particle1_pos = particle1.getPosition();
                const std::array<double, 3>& particle2_pos = particle2.getPosition();
                double charge1 = particle1.getCharge();
                double charge2 = particle2.getCharge();

                double dx = particle2_pos[0] - particle1_pos[0];
                double dy = particle2_pos[1] - particle1_pos[1];
                double dz = particle2_pos[2] - particle1_pos[2];

                // Calculate distance between atoms
                double distance_squared = dx*dx + dy*dy + dz*dz;
                double distance = std::sqrt(distance_squared);

                // Skip if atoms are too close to avoid division by zero
                if (distance < 1e-10) continue;

                // Calculate the magnitude of the electrostatic force using Coulomb's law
                double force_magnitude = coulombConstant * (charge1 * charge2) / distance_squared;

                // Calculate force components by scaling displacement vector
                double force_scale = force_magnitude / distance;
                double force_x = force_scale * dx;
                double force_y = force_scale * dy;
                double force_z = force_scale * dz;

                thread_local_forces_x[particle1_index] -= force_x;
                thread_local_forces_y[particle1_index] -= force_y;
                thread_local_forces_z[particle1_index] -= force_z;

                thread_local_forces_x[particle2_index] += force_x;
                thread_local_forces_y[particle2_index] += force_y;
                thread_local_forces_z[particle2_index] += force_z;

                VERBOSE("forces: %f %f %f\n", force_x, force_y, force_z);
                VERBOSE("positions 1: %f %f %f\n", particle1_pos[0], particle1_pos[1], particle1_pos[2]);
                VERBOSE("positions 2: %f %f %f\n", particle2_pos[0], particle2_pos[1], particle2_pos[2]);

            }
        }

        # pragma omp critical
        for (size_t particle_index; particle_index < num_particles; ++particle_index) {
            VERBOSE("thread_id: %d\n", thread_id);
            VERBOSE("particle_index: %ld\n", particle_index);
            VERBOSE("local forces: %f %f %f\n", thread_local_forces_x[particle_index], thread_local_forces_y[particle_index], thread_local_forces_z[particle_index]);
            VERBOSE("forces: %f %f %f\n", forces[particle_index][0], forces[particle_index][1], forces[particle_index][2]);
            forces[particle_index][0] += thread_local_forces_x[particle_index];
            forces[particle_index][1] += thread_local_forces_y[particle_index];
            forces[particle_index][2] += thread_local_forces_z[particle_index];
        }

        chronometer.end("electrostatic_force", thread_id);
        chronometer.printTiming("electrostatic_force", "ms", thread_id);
    }

    // Combine thread-local forces into the shared forces array
    // #pragma omp parallel for num_threads(num_threads)


    chronometer.end("electrostatic_force_glob");
    chronometer.printTiming("electrostatic_force_glob", "ms");
}

void ElectrostaticForce::compute(System& system,
                                  std::vector<std::array<double, 3>>& forces,
                                  ThreadManager& thread_manager,
                                  Chronometer& chronometer) const {
    VERBOSE("Computing electrostatic forces\n");
    printf("Computing electrostatic forces sequentially\n");

    chronometer.start("electrostatic_force_glob");

    // Create thread-local forces arrays
    size_t num_bonds = bonds.size();
    size_t num_particles = forces.size();

    // std::vector<std::vector<std::array<double, 3>>> local_forces(num_threads,
    //     std::vector<std::array<double, 3>>(forces.size(), {0.0, 0.0, 0.0}));

    // #pragma omp parallel num_threads(num_threads)
    // {
        // int thread_id = omp_get_thread_num();
        chronometer.start("electrostatic_force", 0);
        // std::vector<double> thread_local_forces_x(num_particles, 0.0);
        // std::vector<double> thread_local_forces_y(num_particles, 0.0);
        // std::vector<double> thread_local_forces_z(num_particles, 0.0);
        // std::vector<std::array<double, 3>> thread_local_forces(num_particles, {0.0, 0.0, 0.0});
        
        // #pragma omp for collapse(2)
        for (size_t particle1_index = 0; particle1_index < num_particles; ++particle1_index) {
            for (size_t particle2_index = particle1_index + 1; particle2_index < num_particles; ++particle2_index) {
                Particle& particle1 = system.getParticle(particle1_index);
                Particle& particle2 = system.getParticle(particle2_index);

                const std::array<double, 3>& particle1_pos = particle1.getPosition();
                const std::array<double, 3>& particle2_pos = particle2.getPosition();
                double charge1 = particle1.getCharge();
                double charge2 = particle2.getCharge();

                double dx = particle2_pos[0] - particle1_pos[0];
                double dy = particle2_pos[1] - particle1_pos[1];
                double dz = particle2_pos[2] - particle1_pos[2];

                // Calculate distance between atoms
                double distance_squared = dx*dx + dy*dy + dz*dz;
                double distance = std::sqrt(distance_squared);

                // Skip if atoms are too close to avoid division by zero
                if (distance < 1e-10) continue;

                // Calculate the magnitude of the electrostatic force using Coulomb's law
                double force_magnitude = coulombConstant * (charge1 * charge2) / distance_squared;

                // Calculate force components by scaling displacement vector
                double force_scale = force_magnitude / distance;
                double force_x = force_scale * dx;
                double force_y = force_scale * dy;
                double force_z = force_scale * dz;

                forces[particle1_index][0] -= force_x;
                forces[particle1_index][1] -= force_y;
                forces[particle1_index][2] -= force_z;

                forces[particle2_index][0] += force_x;
                forces[particle2_index][1] += force_y;
                forces[particle2_index][2] += force_z;

                VERBOSE("forces: %f %f %f\n", force_x, force_y, force_z);
                VERBOSE("positions 1: %f %f %f\n", particle1_pos[0], particle1_pos[1], particle1_pos[2]);
                VERBOSE("positions 2: %f %f %f\n", particle2_pos[0], particle2_pos[1], particle2_pos[2]);

            }
        }

        // // # pragma omp critical
        // for (size_t particle_index; particle_index < num_particles; ++particle_index) {
        //     VERBOSE("thread_id: %d\n", 0);
        //     VERBOSE("particle_index: %ld\n", particle_index);
        //     // VERBOSE("local forces: %f %f %f\n", forces[particle_index][0], forces[particle_index][1], forces[particle_index][2]);
        //     VERBOSE("forces: %f %f %f\n", forces[particle_index][0], forces[particle_index][1], forces[particle_index][2]);
        //     forces[particle_index][0] += thread_local_forces_x[particle_index];
        //     forces[particle_index][1] += thread_local_forces_y[particle_index];
        //     forces[particle_index][2] += thread_local_forces_z[particle_index];
        // }

        chronometer.end("electrostatic_force", 0);
        chronometer.printTiming("electrostatic_force", "ms", 0);
    // }

    // Combine thread-local forces into the shared forces array
    // #pragma omp parallel for num_threads(num_threads)


    chronometer.end("electrostatic_force_glob");
    chronometer.printTiming("electrostatic_force_glob", "ms");
}