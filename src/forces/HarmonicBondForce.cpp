#include "HarmonicBondForce.h"
#include <cmath> // for sqrt and pow
// #include <omp.h>    // for OpenMP
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


    for (size_t bondIdx = 0; bondIdx < num_bonds; ++bondIdx) {

        double z = 2.0;
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
}
