#include "HarmonicBondForce.h"
#include <cmath>
#include <omp.h>
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

    printf("Computing harmonic bond forces in parallel\n");
    int num_threads = thread_manager.harmonic_bond_threads;
    size_t num_bonds = bonds.size();
    size_t num_particles = forces.size();
    #pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();
        std::vector<double> thread_local_forces_x(num_particles, 0.0);
        std::vector<double> thread_local_forces_y(num_particles, 0.0);
        std::vector<double> thread_local_forces_z(num_particles, 0.0);

        #pragma omp for
        for (size_t bondIdx = 0; bondIdx < num_bonds; ++bondIdx) {
            const auto& bond = bonds[bondIdx];
            size_t particle1_index = bond.first;
            size_t particle2_index = bond.second;

            Particle& particle1 = system.getParticle(particle1_index);
            Particle& particle2 = system.getParticle(particle2_index);

            const std::array<double, 3>& particle1_pos = particle1.getPosition();
            const std::array<double, 3>& particle2_pos = particle2.getPosition();

            VERBOSE("Bond: %ld %ld\n", particle1_index, particle2_index);

            double dx = particle2_pos[0] - particle1_pos[0];
            double dy = particle2_pos[1] - particle1_pos[1];
            double dz = particle2_pos[2] - particle1_pos[2];
            
            // Calculate distance between atoms
            double distance_squared = dx*dx + dy*dy + dz*dz;
            double distance = std::sqrt(distance_squared);

            // Skip if atoms are too close to avoid division by zero
            if (distance < 1e-10) continue;

            // Calculate the magnitude of the bond force using Coulomb's law
            double displacement = distance - equilibriumDistance;
            double force_magnitude = springConstant * displacement;

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
        }

        # pragma omp critical
        for (size_t particle_index; particle_index < num_particles; ++particle_index) {
            forces[particle_index][0] += thread_local_forces_x[particle_index];
            forces[particle_index][1] += thread_local_forces_y[particle_index];
            forces[particle_index][2] += thread_local_forces_z[particle_index];
        }
    }
}

void HarmonicBondForce::compute(System& system, 
                                    std::vector<std::array<double, 3>>& forces,
                                    ThreadManager& thread_manager,
                                    Chronometer& chronometer) const {

    printf("Computing harmonic bond forces sequentially\n");
    int num_threads = thread_manager.harmonic_bond_threads;
    size_t num_bonds = bonds.size();
    size_t num_particles = forces.size();
    for (size_t bondIdx = 0; bondIdx < num_bonds; ++bondIdx) {
        const auto& bond = bonds[bondIdx];
        size_t particle1_index = bond.first;
        size_t particle2_index = bond.second;

        Particle& particle1 = system.getParticle(particle1_index);
        Particle& particle2 = system.getParticle(particle2_index);

        const std::array<double, 3>& particle1_pos = particle1.getPosition();
        const std::array<double, 3>& particle2_pos = particle2.getPosition();

        VERBOSE("Bond: %ld %ld\n", particle1_index, particle2_index);

        double dx = particle2_pos[0] - particle1_pos[0];
        double dy = particle2_pos[1] - particle1_pos[1];
        double dz = particle2_pos[2] - particle1_pos[2];

        // Calculate distance between atoms
        double distance_squared = dx*dx + dy*dy + dz*dz;
        double distance = std::sqrt(distance_squared);

        // Skip if atoms are too close to avoid division by zero
        if (distance < 1e-10) continue;

        // Calculate the magnitude of the electrostatic force using Coulomb's law
        double displacement = distance - equilibriumDistance;
        double force_magnitude = springConstant * displacement;

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
    }
}