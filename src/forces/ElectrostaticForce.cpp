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

void ElectrostaticForce::compute(System& system,
                                  std::vector<std::array<double, 3>>& forces,
                                  ThreadManager& thread_manager) const {
    auto start_time = std::chrono::high_resolution_clock::now();

    // Create thread-local forces arrays
    int num_threads = thread_manager.electrostatic_threads;
    size_t num_bonds = bonds.size();
    size_t num_particles = forces.size();

    std::vector<std::vector<std::array<double, 3>>> local_forces(num_threads,
        std::vector<std::array<double, 3>>(forces.size(), {0.0, 0.0, 0.0}));

    #pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();
        std::vector<double> thread_local_forces_x(num_particles, 0.0);
        std::vector<double> thread_local_forces_y(num_particles, 0.0);
        std::vector<double> thread_local_forces_z(num_particles, 0.0);
        // std::vector<std::array<double, 3>> thread_local_forces(num_particles, {0.0, 0.0, 0.0});
        
        #pragma omp for
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

                // Accumulate forces in thread-local storage
                // Force on particle1
                // thread_local_forces[particle1_index][0] += force_x;
                // thread_local_forces[particle1_index][1] += force_y;
                // thread_local_forces[particle1_index][2] += force_z;

                // // Equal and opposite force on particle2
                // thread_local_forces[particle2_index][0] -= force_x;
                // thread_local_forces[particle2_index][1] -= force_y;
                // thread_local_forces[particle2_index][2] -= force_z;

                thread_local_forces_x[particle1_index] += force_x;
                thread_local_forces_y[particle1_index] += force_y;
                thread_local_forces_z[particle1_index] += force_z;

                thread_local_forces_x[particle2_index] -= force_x;
                thread_local_forces_y[particle2_index] -= force_y;
                thread_local_forces_z[particle2_index] -= force_z;

                // forces[particle1_index][0] += force_x;
                // forces[particle1_index][1] += force_y;
                // forces[particle1_index][2] += force_z;

                // // Equal and opposite force on particle2
                // forces[particle2_index][0] -= force_x;
                // forces[particle2_index][1] -= force_y;
                // forces[particle2_index][2] -= force_z;
                VERBOSE("forces: %f %f %f\n", force_x, force_y, force_z);
                VERBOSE("positions 1: %f %f %f\n", particle1_pos[0], particle1_pos[1], particle1_pos[2]);
                VERBOSE("positions 2: %f %f %f\n", particle2_pos[0], particle2_pos[1], particle2_pos[2]);

            }
        }
        
        # pragma omp critical
        for (size_t thread_index, particle_index = 0; particle_index < num_particles; ++particle_index) {
            forces[particle_index][0] += thread_local_forces_x[particle_index];
            forces[particle_index][1] += thread_local_forces_y[particle_index];
            forces[particle_index][2] += thread_local_forces_z[particle_index];
        }




            // for (size_t i = 0; i < forces.size(); ++i) {
            //     for (size_t k = 0; k < 3; ++k) {
            //         for (int t = 0; t < num_threads; ++t) {
            //             printf("local forces: %f %f %f\n", local_forces[t][i][0], local_forces[t][i][1], local_forces[t][i][2]);
            //             forces[i][k] += thread_local_forces[t][i][k];
            //         }
            //     }
            // }
    }

    // Combine thread-local forces into the shared forces array
    // #pragma omp parallel for num_threads(num_threads)


    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    printf("Time taken by electrostatic force computation: %d ms\n", duration.count());
}