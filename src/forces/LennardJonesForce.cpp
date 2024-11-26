#include "LennardJonesForce.h"
#include "../logging/Verbose.h"
#include <omp.h>

LennardJonesForce::LennardJonesForce(double epsilon, double sigma)
    : epsilon(epsilon), sigma(sigma) {}

void LennardJonesForce::compute_parallel(System& system, std::vector<std::array<double, 3>>& forces,
                                    ThreadManager& thread_manager,
                                    Chronometer& chronometer) const {
    printf("Computing dispersion forces in parallel\n");
    int num_threads = thread_manager.lennard_jones_threads;
    size_t num_particles = forces.size();
    #pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();
        std::vector<double> thread_local_forces_x(num_particles, 0.0);
        std::vector<double> thread_local_forces_y(num_particles, 0.0);
        std::vector<double> thread_local_forces_z(num_particles, 0.0);
        
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

                VERBOSE("Dispersion interaction: %ld %ld\n", particle1_index, particle2_index);

                // Calculate distance between atoms
                double distance_squared = dx*dx + dy*dy + dz*dz;
                double distance = std::sqrt(distance_squared);

                // Skip if atoms are too close to avoid division by zero
                if (distance < 1e-10) continue;

                // Calculate Lennard-Jones 6-12 potential terms
                double sr2 = sigma * sigma / distance_squared; // (σ / r)^2
                double sr6 = sr2 * sr2 * sr2;                 // (σ / r)^6
                double sr12 = sr6 * sr6;                      // (σ / r)^12

                // Calculate the magnitude of the Lennard-Jones force
                double force_magnitude = 24 * epsilon * (2 * sr12 - sr6) / distance_squared;

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
        }

        # pragma omp critical
        for (size_t particle_index; particle_index < num_particles; ++particle_index) {
            forces[particle_index][0] += thread_local_forces_x[particle_index];
            forces[particle_index][1] += thread_local_forces_y[particle_index];
            forces[particle_index][2] += thread_local_forces_z[particle_index];
        }
    }
}


void LennardJonesForce::compute(System& system,
                                  std::vector<std::array<double, 3>>& forces,
                                  ThreadManager& thread_manager,
                                  Chronometer& chronometer) const {
    printf("Computing dispersion forces sequentially\n");
    size_t num_particles = forces.size();

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

            VERBOSE("Electrostatic nteraction: %ld %ld\n", particle1_index, particle2_index);

            // Calculate distance between atoms
            double distance_squared = dx*dx + dy*dy + dz*dz;
            double distance = std::sqrt(distance_squared);

            // Skip if atoms are too close to avoid division by zero
            if (distance < 1e-10) continue;

            // Calculate the magnitude of the electrostatic force using Coulomb's law
            double sr2 = sigma * sigma / distance_squared; // (σ / r)^2
            double sr6 = sr2 * sr2 * sr2;                 // (σ / r)^6
            double sr12 = sr6 * sr6;                      // (σ / r)^12

            // Calculate the magnitude of the Lennard-Jones force
            double force_magnitude = 24 * epsilon * (2 * sr12 - sr6) / distance_squared;

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
}