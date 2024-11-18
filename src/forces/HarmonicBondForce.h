#ifndef HARMONIC_BOND_FORCE_H
#define HARMONIC_BOND_FORCE_H

#include "Force.h"
// #include "System.h"
#include <utility>
#include <vector>

/**
 * Represents a harmonic bond force between pairs of particles.
 */
class HarmonicBondForce : public Force {
public:
    HarmonicBondForce(double springConstant, double equilibriumDistance);

    /**
     * Adds a bond between two particles by their indices.
     * @param particle1Index - Index of the first particle.
     * @param particle2Index - Index of the second particle.
     */
    void addBond(size_t particle1Index, size_t particle2Index);

    /**
     * Compute the harmonic bond forces and update the force array.
     * @param system - The system containing particles.
     * @param forces - A reference to the forces array to update.
     */
    void compute(System& system, std::vector<std::array<double, 3>>& forces) const override;

private:
    double springConstant;
    double equilibriumDistance;
    std::vector<std::pair<size_t, size_t>> bonds; // List of particle pairs for bonds
};

#endif // HARMONIC_BOND_FORCE_H
