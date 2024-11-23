#include "LennardJonesForce.h"

LennardJonesForce::LennardJonesForce(double epsilon, double sigma)
    : epsilon(epsilon), sigma(sigma) {}

void LennardJonesForce::compute(System& system, std::vector<std::array<double, 3>>& forces,
                                    ThreadManager& tread_manager,
                                    Chronometer& chronometer) const {
    size_t numParticles = system.getNumParticles();

    for (size_t i = 0; i < numParticles; ++i) {
        Particle& particle1 = system.getParticle(i);
        const std::array<double, 3>& pos1 = particle1.getPosition();

        for (size_t j = i + 1; j < numParticles; ++j) {
            Particle& particle2 = system.getParticle(j);
            const std::array<double, 3>& pos2 = particle2.getPosition();

            std::array<double, 3> r = {0.0, 0.0, 0.0};
            double rSquared = 0.0;
            for (int k = 0; k < 3; ++k) {
                r[k] = pos2[k] - pos1[k];
                rSquared += r[k] * r[k];
            }
            double rMagnitude = std::sqrt(rSquared);
            if (rMagnitude == 0.0) continue;

            double sr6 = std::pow(sigma / rMagnitude, 6);
            double forceMagnitude = 24 * epsilon * (2 * sr6 * sr6 - sr6) / rSquared;

            for (int k = 0; k < 3; ++k) {
                double force = forceMagnitude * r[k];
                forces[i][k] += force;
                forces[j][k] -= force;
            }
        }
    }
}
