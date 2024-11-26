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

void HarmonicBondForce::compute(System& system, std::vector<std::array<double, 3>>& forces) const {
    for (const auto& bond : bonds) {
        size_t i = bond.first;
        size_t j = bond.second;

        Particle& particle1 = system.getParticle(i);
        Particle& particle2 = system.getParticle(j);

        const std::array<double, 3>& pos1 = particle1.getPosition();
        const std::array<double, 3>& pos2 = particle2.getPosition();

        // Calculate the distance vector and its magnitude
        std::array<double, 3> r = {0.0, 0.0, 0.0};
        double rSquared = 0.0;

        // Calculate displacement vector and squared distance
        for (int k = 0; k < 3; ++k) {
            r[k] = pos2[k] - pos1[k];
            rSquared += r[k] * r[k];
        }

        // Calculate actual distance between particles
        double rMagnitude = std::sqrt(rSquared);

        // Skip if particles are at the same position to avoid division by zero
        if (rMagnitude < 1e-10) continue;

        // Calculate the force using the harmonic potential formula:
        // F = -k(r - r₀)r̂
        // where k is the spring constant, r is the current distance,
        // r₀ is the equilibrium distance, and r̂ is the unit vector
        double deltaR = rMagnitude - equilibriumDistance;
        double forceMagnitude = springConstant * deltaR;

        // Convert to force per component by multiplying by the normalized direction vector
        for (int k = 0; k < 3; ++k) {
            // Divide by rMagnitude to normalize the direction vector
            double force = forceMagnitude * (r[k] / rMagnitude);

            // Apply equal and opposite forces to both particles (Newton's Third Law)
            forces[i][k] += force;  // Force on particle i
            forces[j][k] -= force;  // Equal and opposite force on particle j
        }
    }
}
