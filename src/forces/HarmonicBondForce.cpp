#include "HarmonicBondForce.h"
#include <cmath> // for sqrt and pow
#include <omp.h>    // for OpenMP
#include "../logging/Verbose.h"
#include <unistd.h>
#include <iostream>
#include <chrono>
#include "../thread_manager/ThreadManager.h"

HarmonicBondForce::HarmonicBondForce(double springConstant, double equilibriumDistance)
    : springConstant(springConstant), equilibriumDistance(equilibriumDistance) {}

void HarmonicBondForce::addBond(size_t particle1Index, size_t particle2Index) {
    bonds.emplace_back(particle1Index, particle2Index);
}

void HarmonicBondForce::compute_parallel(System& system, 
                                    std::vector<std::array<double, 3>>& forces,
                                    ThreadManager& thread_manager,
                                    Chronometer& chronometer) const {
    // extern ThreadConfig thread_config;
    // printf("harmonic_bond_threads: %d\n", thread_manager.harmonic_bond_threads);

    // auto start_time = std::chrono::high_resolution_clock::now();

    // Create thread-local forces arrays
    int num_threads = thread_manager.harmonic_bond_threads;
    size_t num_bonds = bonds.size();
    // printf("num_bonds: %ld\n", num_bonds);
    size_t num_particles = forces.size();
    // printf("num_particles: %ld\n", num_particles);


    std::vector<std::vector<std::array<double, 3>>> local_forces(num_threads, 
    std::vector<std::array<double, 3>>(forces.size(), {0.0, 0.0, 0.0}));


    
    #pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();
        printf("thread_id: %d\n", thread_id);
        // chronometer.start("electrostatic_force", thread_id);
        std::vector<double> thread_local_forces_x(num_particles, 0.0);
        std::vector<double> thread_local_forces_y(num_particles, 0.0);
        std::vector<double> thread_local_forces_z(num_particles, 0.0);

        #pragma omp for
        for (size_t bondIdx = 0; bondIdx < num_bonds; ++bondIdx) {
            // sleep 1 sec with c style code

            // #include <unistd.h>
            double z = 2.0;
            // for (int i = 0; i < 1000; i++) {
            //     z += std::sqrt(z+1);
            // }
            // printf("bond %d Thread %d: %d\n", bondIdx, 0, z);
            const auto& bond = bonds[bondIdx];
            size_t particle1_index = bond.first;
            size_t particle2_index = bond.second;

            Particle& particle1 = system.getParticle(particle1_index);
            Particle& particle2 = system.getParticle(particle2_index);

            const std::array<double, 3>& particle1_pos = particle1.getPosition();
            const std::array<double, 3>& particle2_pos = particle2.getPosition();

            double dx = particle2_pos[0] - particle1_pos[0];
            double dy = particle2_pos[1] - particle1_pos[1];
            double dz = particle2_pos[2] - particle1_pos[2];

            // Calculate distance between atoms
            double distance_squared = dx*dx + dy*dy + dz*dz;
            double distance = std::sqrt(distance_squared);

            // Skip if atoms are too close to avoid division by zero
            if (distance < 1e-10) continue;

            // Calculate magnitude of harmonic force
            double displacement = distance - equilibriumDistance;
            double force_magnitude = springConstant * displacement;

            // Calculate force components by scaling displacement vector
            double force_scale = force_magnitude / distance;
            double force_x = force_scale * dx;
            double force_y = force_scale * dy;
            double force_z = force_scale * dz;

            // Accumulate forces in thread-local storage
            // Force on atom1
            thread_local_forces_x[particle1_index] -= force_x;
            thread_local_forces_y[particle1_index] -= force_y;
            thread_local_forces_z[particle1_index] -= force_z;

            thread_local_forces_x[particle2_index] += force_x;
            thread_local_forces_y[particle2_index] += force_y;
            thread_local_forces_z[particle2_index] += force_z;

        }

        // auto end_loop_time = std::chrono::high_resolution_clock::now();
        // auto loop_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_loop_time - start_loop_time);
        // if (thread_id == 0) {
            // std::cout << "Time taken by loop iterations in parallel region: " 
            //           << loop_duration.count() << " ms" << std::endl;
            // printf("time taken by loop iterations in parallel region: %d ms\n", loop_duration.count());
        // }
        # pragma omp critical
        for (size_t particle_index; particle_index < num_particles; ++particle_index) {
            // VERBOSE("thread_id: %d\n", thread_id);
            VERBOSE("particle_index: %ld\n", particle_index);
            VERBOSE("local forces: %f %f %f\n", thread_local_forces_x[particle_index], thread_local_forces_y[particle_index], thread_local_forces_z[particle_index]);
            VERBOSE("forces: %f %f %f\n", forces[particle_index][0], forces[particle_index][1], forces[particle_index][2]);
            forces[particle_index][0] += thread_local_forces_x[particle_index];
            forces[particle_index][1] += thread_local_forces_y[particle_index];
            forces[particle_index][2] += thread_local_forces_z[particle_index];
        }
    }

    // // Combine thread-local forces into the shared forces array
    // // #pragma omp parallel for num_threads(num_threads)
    // for (size_t i = 0; i < forces.size(); ++i) {
    //     for (size_t k = 0; k < 3; ++k) {
    //         for (int t = 0; t < num_threads; ++t) {
    //             forces[i][k] += local_forces[t][i][k];
    //         }
    //     }
    // }

    // auto end_time = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    // printf("Time taken by harmonic bond force computation: %d ms\n", duration.count());
}

void HarmonicBondForce::compute(System& system, 
                                            std::vector<std::array<double, 3>>& forces,
                                            ThreadManager& thread_manager,
                                            Chronometer& chronometer) const {

    int num_threads = thread_manager.harmonic_bond_threads;
    size_t num_bonds = bonds.size();
    // printf("num_bonds: %ld\n", num_bonds);
    size_t num_particles = forces.size();
    // printf("num_particles: %ld\n", num_particles);


    std::vector<std::vector<std::array<double, 3>>> local_forces(num_threads, 
    std::vector<std::array<double, 3>>(forces.size(), {0.0, 0.0, 0.0}));


    
    // #pragma omp parallel num_threads(num_threads)
    // {
    //     int thread_id = omp_get_thread_num();
    //     printf("thread_id: %d\n", thread_id);
    //     // chronometer.start("electrostatic_force", thread_id);
    //     std::vector<double> thread_local_forces_x(num_particles, 0.0);
    //     std::vector<double> thread_local_forces_y(num_particles, 0.0);
    //     std::vector<double> thread_local_forces_z(num_particles, 0.0);

    //     #pragma omp for
    for (size_t bondIdx = 0; bondIdx < num_bonds; ++bondIdx) {
        // sleep 1 sec with c style code

        // #include <unistd.h>
        double z = 2.0;
        // for (int i = 0; i < 1000; i++) {
        //     z += std::sqrt(z+1);
        // }
        // printf("bond %d Thread %d: %d\n", bondIdx, 0, z);
        const auto& bond = bonds[bondIdx];
        size_t particle1_index = bond.first;
        size_t particle2_index = bond.second;

        Particle& particle1 = system.getParticle(particle1_index);
        Particle& particle2 = system.getParticle(particle2_index);

        const std::array<double, 3>& particle1_pos = particle1.getPosition();
        const std::array<double, 3>& particle2_pos = particle2.getPosition();

        double dx = particle2_pos[0] - particle1_pos[0];
        double dy = particle2_pos[1] - particle1_pos[1];
        double dz = particle2_pos[2] - particle1_pos[2];

        // Calculate distance between atoms
        double distance_squared = dx*dx + dy*dy + dz*dz;
        double distance = std::sqrt(distance_squared);

        // Skip if atoms are too close to avoid division by zero
        if (distance < 1e-10) continue;

        // Calculate magnitude of harmonic force
        double displacement = distance - equilibriumDistance;
        double force_magnitude = springConstant * displacement;

        // Calculate force components by scaling displacement vector
        double force_scale = force_magnitude / distance;
        double force_x = force_scale * dx;
        double force_y = force_scale * dy;
        double force_z = force_scale * dz;

        // Accumulate forces in thread-local storage
        // Force on atom1
        forces[particle1_index][0] -= force_x;
        forces[particle1_index][1] -= force_y;
        forces[particle1_index][2] -= force_z;

        forces[particle2_index][0] += force_x;
        forces[particle2_index][1] += force_y;
        forces[particle2_index][2] += force_z;

    }

        // // auto end_loop_time = std::chrono::high_resolution_clock::now();
        // // auto loop_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_loop_time - start_loop_time);
        // // if (thread_id == 0) {
        //     // std::cout << "Time taken by loop iterations in parallel region: " 
        //     //           << loop_duration.count() << " ms" << std::endl;
        //     // printf("time taken by loop iterations in parallel region: %d ms\n", loop_duration.count());
        // // }
        // # pragma omp critical
        // for (size_t particle_index; particle_index < num_particles; ++particle_index) {
        //     VERBOSE("thread_id: %d\n", thread_id);
        //     VERBOSE("particle_index: %ld\n", particle_index);
        //     VERBOSE("local forces: %f %f %f\n", thread_local_forces_x[particle_index], thread_local_forces_y[particle_index], thread_local_forces_z[particle_index]);
        //     VERBOSE("forces: %f %f %f\n", forces[particle_index][0], forces[particle_index][1], forces[particle_index][2]);
        //     forces[particle_index][0] += thread_local_forces_x[particle_index];
        //     forces[particle_index][1] += thread_local_forces_y[particle_index];
        //     forces[particle_index][2] += thread_local_forces_z[particle_index];
        // }
}

    // // Combine thread-local forces into the shared forces array
    // // #pragma omp parallel for num_threads(num_threads)
    // for (size_t i = 0; i < forces.size(); ++i) {
    //     for (size_t k = 0; k < 3; ++k) {
    //         for (int t = 0; t < num_threads; ++t) {
    //             forces[i][k] += local_forces[t][i][k];
    //         }
    //     }
    // }

    // auto end_time = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    // printf("Time taken by harmonic bond force computation: %d ms\n", duration.count());

