#include "HarmonicBondForce.h"
#include <cmath> // for sqrt and pow
#include <omp.h>    // for OpenMP
#include "../logging/Verbose.h"

HarmonicBondForce::HarmonicBondForce(double springConstant, double equilibriumDistance)
    : springConstant(springConstant), equilibriumDistance(equilibriumDistance) {}

void HarmonicBondForce::addBond(size_t particle1Index, size_t particle2Index) {
    bonds.emplace_back(particle1Index, particle2Index);
}

// void HarmonicBondForce::compute(System& system, std::vector<std::array<double, 3>>& forces) const {
//     // std::vector<std::array<double, 3>> local_forces(forces.size(), {0.0, 0.0, 0.0}); // Create a local forces array
//     #pragma omp parallel for num_threads(4)
//     for (const auto& bond : bonds) {
//         size_t i = bond.first;
//         size_t j = bond.second;

//         Particle& particle1 = system.getParticle(i);
//         Particle& particle2 = system.getParticle(j);

//         const std::array<double, 3>& pos1 = particle1.getPosition();
//         const std::array<double, 3>& pos2 = particle2.getPosition();

//         // Calculate the distance vector and its magnitude
//         std::array<double, 3> r = {0.0, 0.0, 0.0};
//         double rSquared = 0.0;
        
//         // Calculate displacement vector and squared distance
//         for (int k = 0; k < 3; ++k) {
//             r[k] = pos2[k] - pos1[k];
//             rSquared += r[k] * r[k];
//         }
        
//         // Calculate actual distance between particles
//         double rMagnitude = std::sqrt(rSquared);

//         // Skip if particles are at the same position to avoid division by zero
//         if (rMagnitude < 1e-10) continue;

//         // Calculate the force using the harmonic potential formula:
//         // F = -k(r - r₀)r̂
//         // where k is the spring constant, r is the current distance,
//         // r₀ is the equilibrium distance, and r̂ is the unit vector
//         double deltaR = rMagnitude - equilibriumDistance;
//         double forceMagnitude = springConstant * deltaR;
        
//         // Convert to force per component by multiplying by the normalized direction vector
//         for (int k = 0; k < 3; ++k) {
//             // Divide by rMagnitude to normalize the direction vector
//             double force = forceMagnitude * (r[k] / rMagnitude);
            
//             // Apply equal and opposite forces to both particles (Newton's Third Law)
//             // #pragma omp atomic
//             forces[i][k] += force;  // Force on particle i
//             // #pragma omp atomic
//             forces[j][k] -= force;  // Equal and opposite force on particle j
//         }
//     }
// }

void HarmonicBondForce::compute(System& system, std::vector<std::array<double, 3>>& forces) const {
    // Create thread-local forces arrays
    int num_threads = 1;
    std::vector<std::vector<std::array<double, 3>>> local_forces(num_threads, 
        std::vector<std::array<double, 3>>(forces.size(), {0.0, 0.0, 0.0}));
    
    #pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num(); // Get the thread ID for the current thread
        auto& thread_local_forces = local_forces[thread_id]; // Reference to this thread's local forces

        #pragma omp for
        for (size_t bondIdx = 0; bondIdx < bonds.size(); ++bondIdx) {
            const auto& bond = bonds[bondIdx];
            size_t particle1_idx = bond.first;
            size_t particle2_idx = bond.second;

            Particle& particle1 = system.getParticle(particle1_idx);
            Particle& particle2 = system.getParticle(particle2_idx);

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
            thread_local_forces[particle1_idx][0] += force_x;
            thread_local_forces[particle1_idx][1] += force_y;
            thread_local_forces[particle1_idx][2] += force_z;

            // Equal and opposite force on particle2
            thread_local_forces[particle2_idx][0] -= force_x;
            thread_local_forces[particle2_idx][1] -= force_y;
            thread_local_forces[particle2_idx][2] -= force_z;
            }
    }

    // Combine thread-local forces into the shared forces array
    for (int t = 0; t < num_threads; ++t) {
        for (size_t i = 0; i < forces.size(); ++i) {
            for (int k = 0; k < 3; ++k) {
                forces[i][k] += local_forces[t][i][k];
            }
        }
    }
}

